// A spatial hash grid for rigid bodies with angular trajectories.
#include "rigid_body_hash_grid.hpp"

#include <cmath>
#include <iostream>

#include <tbb/parallel_invoke.h>

#include <interval/interval.hpp>
#include <logger.hpp>

#include <io/serialize_json.hpp>
#include <nlohmann/json.hpp>

namespace ccd {

/// Resize to fit a static scene
void RigidBodyHashGrid::resize(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const std::vector<int>& body_ids,
    const double inflation_radius)
{
    int dim = bodies.dim();
    Eigen::VectorX3d min, max;
    min.setConstant(dim, +std::numeric_limits<double>::infinity());
    max.setConstant(dim, -std::numeric_limits<double>::infinity());

    for (int i : body_ids) {
        Eigen::MatrixXd V = bodies[i].world_vertices(poses[i]);
        min = min.cwiseMin(V.colwise().minCoeff().transpose());
        max = max.cwiseMax(V.colwise().maxCoeff().transpose());
    }

    min.array() -= inflation_radius;
    min.array() += inflation_radius;

    // TODO: this may not be well scaled depending on the body_ids
    double cell_size = bodies.average_edge_length;

    HashGrid::resize(min, max, cell_size);
}

void compute_scene_conservative_bbox(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const std::vector<int>& body_ids,
    Eigen::VectorX3d& min,
    Eigen::VectorX3d& max)
{
    int dim = bodies.dim();
    min.setConstant(dim, +std::numeric_limits<double>::infinity());
    max.setConstant(dim, -std::numeric_limits<double>::infinity());

    for (int i : body_ids) {
        const auto& body = bodies[i];

        const Eigen::VectorX3d& p0 = poses_t0[i].position;
        const Eigen::VectorX3d& p1 = poses_t1[i].position;
        const Eigen::VectorX3d& r0 = poses_t0[i].rotation;
        const Eigen::VectorX3d& r1 = poses_t1[i].rotation;

        // If the body is not rotating then just use the linearized trajectory
        if ((r0.array() == r1.array()).all()) {
            Eigen::MatrixXd V0 = body.world_vertices(poses_t0[i]);
            min = min.cwiseMin(V0.colwise().minCoeff().transpose());
            max = max.cwiseMax(V0.colwise().maxCoeff().transpose());
            Eigen::MatrixXd V1 = body.world_vertices(poses_t1[i]);
            min = min.cwiseMin(V1.colwise().minCoeff().transpose());
            max = max.cwiseMax(V1.colwise().maxCoeff().transpose());
        } else {
            // Use the maximum radius of the body to bound all rotations
            min = min.cwiseMin((p0.array() - body.r_max).matrix());
            max = max.cwiseMax((p0.array() + body.r_max).matrix());
            min = min.cwiseMin((p1.array() - body.r_max).matrix());
            max = max.cwiseMax((p1.array() + body.r_max).matrix());
        }
    }

    // Account for rounding controls of interval arithmetic
    double diag_len = (max - min).norm();
    min.array() -= 1e-8 * diag_len;
    max.array() += 1e-8 * diag_len;
}

void RigidBodyHashGrid::resize(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const std::vector<int>& body_ids,
    const double inflation_radius)
{
    Eigen::VectorX3d min(bodies.dim()), max(bodies.dim());
    compute_scene_conservative_bbox(
        bodies, poses_t0, poses_t1, body_ids, min, max);
    min.array() -= inflation_radius;
    min.array() += inflation_radius;

    // Compute the average displacement length
    // NOTE: this is a linearized displacement
    Eigen::MatrixXd V0 = bodies.world_vertices(poses_t0);
    Eigen::MatrixXd V1 = bodies.world_vertices(poses_t1);
    double average_displacement_length = (V1 - V0).rowwise().norm().mean();

    // TODO: this may not be well scaled depending on the body_ids
    double average_edge_length = bodies.average_edge_length;

    double cell_size =
        std::max(average_displacement_length, average_edge_length);

    HashGrid::resize(min, max, cell_size);
}

/// Add dynamic bodies
void RigidBodyHashGrid::addBodies(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const std::vector<int>& body_ids,
    const double inflation_radius)
{
    for (int id : body_ids) {
        Eigen::MatrixXd V = bodies[id].world_vertices(poses[id]);
        tbb::parallel_invoke(
            [&]() {
                long v0i = bodies.m_body_vertex_id[id];
                for (int i = 0; i < V.rows(); i++) {
                    this->addVertex(
                        V.row(i), V.row(i), v0i + i, inflation_radius);
                }
            },
            [&]() {
                const Eigen::MatrixXi& E = bodies[id].edges;
                long e0i = bodies.m_body_edge_id[id];
                for (int i = 0; i < E.rows(); i++) {
                    this->addEdge(
                        V.row(E(i, 0)), V.row(E(i, 1)), //
                        V.row(E(i, 0)), V.row(E(i, 1)), //
                        e0i + i, inflation_radius);
                }
            },
            [&]() {
                const Eigen::MatrixXi& F = bodies[id].faces;
                long f0i = bodies.m_body_face_id[id];
                for (int i = 0; i < F.rows(); i++) {
                    this->addFace(
                        V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)), //
                        V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)), //
                        f0i + i, inflation_radius);
                }
            });
    }
}

void RigidBodyHashGrid::compute_vertices_intervals(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<Interval>& poses_t0,
    const physics::Poses<Interval>& poses_t1,
    const std::vector<int>& body_ids,
    Eigen::MatrixXI& vertices,
    double inflation_radius) const
{
    vertices.setConstant(
        bodies.num_vertices(), bodies.dim(), Interval::empty());
    for (int i : body_ids) {
        Eigen::MatrixXI V;
        int n_subs = compute_vertices_intervals(
            bodies[i], poses_t0[i], poses_t1[i], V, inflation_radius,
            Interval(0, 1), 1);
        if (n_subs) {
            spdlog::trace("nsubs={:d}", n_subs);
        }
        vertices.middleRows(bodies.m_body_vertex_id[i], V.rows()) = V;
    }
}

int RigidBodyHashGrid::compute_vertices_intervals(
    const physics::RigidBody& body,
    const physics::Pose<Interval>& pose_t0,
    const physics::Pose<Interval>& pose_t1,
    Eigen::MatrixX<Interval>& vertices,
    double inflation_radius, // Only used for fit check
    const Interval& t,
    int force_subdivision) const
{
    if (force_subdivision <= 0) {
        physics::Pose<Interval> pose =
            physics::Pose<Interval>::interpolate(pose_t0, pose_t1, t);
        vertices = body.world_vertices(pose);
        // Check that the vertex intervals fit inside the scene bbox
        bool fits = true;
        for (int i = 0; i < vertices.rows(); i++) {
            for (int j = 0; j < vertices.cols(); j++) {
                if (vertices(i, j).lower() - inflation_radius < m_domainMin(j)
                    || vertices(i, j).upper() + inflation_radius
                        > m_domainMax(j)) {
                    fits = false;
                    break;
                }
            }
        }

        if (fits) {
            return 0;
        }
    } else {
        vertices.resizeLike(body.vertices);
    }

    force_subdivision--;

    // If the vertices' intervals are outside the scene bbox, then split t in
    // hopes that a smaller interval will be more accurate.
    std::pair<Interval, Interval> t_halves = bisect(t);
    Eigen::MatrixXI V_first, V_second;
    int n_subs0 = compute_vertices_intervals(
        body, pose_t0, pose_t1, V_first, inflation_radius, t_halves.first,
        force_subdivision);
    int n_subs1 = compute_vertices_intervals(
        body, pose_t0, pose_t1, V_second, inflation_radius, t_halves.second,
        force_subdivision);
    assert(vertices.rows() == V_first.rows());
    assert(vertices.rows() == V_second.rows());
    assert(vertices.cols() == V_first.cols());
    assert(vertices.cols() == V_second.cols());
    // Take the hull of the two halves of the vertices' trajectories
    for (int i = 0; i < vertices.rows(); i++) {
        for (int j = 0; j < vertices.cols(); j++) {
            vertices(i, j) = hull(V_first(i, j), V_second(i, j));
            assert(vertices(i, j).lower() - inflation_radius >= m_domainMin(j));
            assert(vertices(i, j).upper() + inflation_radius <= m_domainMax(j));
        }
    }

    return std::max(n_subs0, n_subs1) + 1;
}

template <typename Derived>
inline ipc::AABB
intervals_to_AABB(const Eigen::MatrixBase<Derived>& x, double inflation_radius)
{
    assert(x.rows() == 1 || x.cols() == 1);
    Eigen::VectorX3d min(x.size());
    Eigen::VectorX3d max(x.size());
    for (int i = 0; i < x.size(); i++) {
        if (empty(x(i))) {
            throw "interval is empty";
        }
        min(i) = x(i).lower() - inflation_radius;
        max(i) = x(i).upper() + inflation_radius;
    }
    return ipc::AABB(min, max);
}

void RigidBodyHashGrid::addBodies(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const std::vector<int>& body_ids,
    const double inflation_radius)
{
    assert(bodies.num_bodies() == poses_t0.size());
    assert(poses_t0.size() == poses_t1.size());

    Eigen::MatrixX<Interval> vertices;
    compute_vertices_intervals(
        bodies, physics::cast<double, Interval>(poses_t0),
        physics::cast<double, Interval>(poses_t1), body_ids, vertices,
        inflation_radius);

    // Create a bounding box for all vertices
    std::vector<ipc::AABB> vertices_aabb;
    vertices_aabb.resize(vertices.rows());
    std::vector<bool> is_vertex_included;
    is_vertex_included.resize(vertices.rows(), true);
    for (long i = 0; i < vertices.rows(); i++) {
        try {
            vertices_aabb[i] =
                intervals_to_AABB(vertices.row(i), inflation_radius);
        } catch (...) {
            is_vertex_included[i] = false;
        }
    }

    auto is_edge_include = [&](size_t ei) {
        return is_vertex_included[bodies.m_edges(ei, 0)]
            && is_vertex_included[bodies.m_edges(ei, 1)];
    };
    auto is_face_include = [&](size_t fi) {
        return is_vertex_included[bodies.m_faces(fi, 0)]
            && is_vertex_included[bodies.m_faces(fi, 1)]
            && is_vertex_included[bodies.m_faces(fi, 2)];
    };

    tbb::parallel_invoke(
        // Add all vertices of the bodies
        [&]() {
            for (long i = 0; i < vertices.rows(); i++) {
                if (is_vertex_included[i]) {
                    this->addElement(vertices_aabb[i], i, this->m_vertexItems);
                }
            }
        },

        // Add all edge of the bodies
        [&]() {
            for (long i = 0; i < bodies.m_edges.rows(); i++) {
                if (is_edge_include(i)) {
                    this->addElement(
                        ipc::AABB(
                            vertices_aabb[bodies.m_edges(i, 0)],
                            vertices_aabb[bodies.m_edges(i, 1)]),
                        i, this->m_edgeItems);
                }
            }
        },

        // Add all faces of the bodies
        [&]() {
            for (long i = 0; i < bodies.m_faces.rows(); i++) {
                if (is_face_include(i)) {
                    this->addElement(
                        ipc::AABB(
                            ipc::AABB(
                                vertices_aabb[bodies.m_faces(i, 0)],
                                vertices_aabb[bodies.m_faces(i, 1)]),
                            vertices_aabb[bodies.m_faces(i, 2)]),
                        i, this->m_faceItems);
                }
            }
        });
}

} // namespace ccd
