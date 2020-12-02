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

void compute_scene_conservative_bbox(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    // const std::vector<int>& body_ids,
    Eigen::VectorX3d& min,
    Eigen::VectorX3d& max)
{
    int dim = bodies.dim();
    min.setConstant(dim, +std::numeric_limits<double>::infinity());
    max.setConstant(dim, -std::numeric_limits<double>::infinity());

    for (int i = 0; i < bodies.m_rbs.size(); i++) {
        const auto& body = bodies.m_rbs[i];

        Eigen::VectorX3d p0 = poses_t0[i].position, p1 = poses_t1[i].position;
        Eigen::VectorX3d r0 = poses_t0[i].rotation, r1 = poses_t1[i].rotation;

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
    compute_scene_conservative_bbox(bodies, poses_t0, poses_t1, min, max);
    min.array() -= inflation_radius;
    min.array() += inflation_radius;

    // Compute the average displacement length
    // NOTE: this is a linearized displacement
    Eigen::MatrixXd V0 = bodies.world_vertices(poses_t0);
    Eigen::MatrixXd V1 = bodies.world_vertices(poses_t1);
    double average_displacement_length = (V1 - V0).rowwise().norm().mean();

    double average_edge_length = bodies.average_edge_length;

    double cell_size =
        std::max(average_displacement_length, average_edge_length);

    HashGrid::resize(min, max, cell_size);
}

void RigidBodyHashGrid::compute_vertices_intervals(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<Interval>& poses_t0,
    const physics::Poses<Interval>& poses_t1,
    // const std::vector<int>& body_ids,
    Eigen::MatrixXI& vertices,
    double inflation_radius) const
{
    vertices.resize(bodies.num_vertices(), bodies.dim());
    int start_i = 0;
    for (int i = 0; i < bodies.num_bodies(); i++) {
        Eigen::MatrixXI V;
        int n_subs = compute_vertices_intervals(
            bodies.m_rbs[i], poses_t0[i], poses_t1[i], V, inflation_radius);
        if (n_subs) {
            spdlog::trace("nsubs={:d}", n_subs);
        }
        vertices.middleRows(start_i, V.rows()) = V;
        start_i += V.rows();
    }
}

int RigidBodyHashGrid::compute_vertices_intervals(
    const physics::RigidBody& body,
    const physics::Pose<Interval>& pose_t0,
    const physics::Pose<Interval>& pose_t1,
    // const std::vector<int>& body_ids,
    Eigen::MatrixX<Interval>& vertices,
    double inflation_radius, // Only used for fit check
    const Interval& t) const
{
    physics::Pose<Interval> pose =
        physics::Pose<Interval>::interpolate(pose_t0, pose_t1, t);
    vertices = body.world_vertices(pose);
    // Check that the vertex intervals fit inside the scene bbox
    bool fits = true;
    for (int i = 0; i < vertices.rows(); i++) {
        for (int j = 0; j < vertices.cols(); j++) {
            if (vertices(i, j).lower() - inflation_radius < m_domainMin(j)
                || vertices(i, j).upper() + inflation_radius > m_domainMax(j)) {
                fits = false;
                break;
            }
        }
    }

    if (fits) {
        return 0;
    }

    // If the vertices' intervals are outside the scene bbox, then split t in
    // hopes that a smaller interval will be more accurate.
    std::pair<Interval, Interval> t_halves = bisect(t);
    Eigen::MatrixXI V_first, V_second;
    int n_subs0 = compute_vertices_intervals(
        body, pose_t0, pose_t1, V_first, inflation_radius, t_halves.first);
    int n_subs1 = compute_vertices_intervals(
        body, pose_t0, pose_t1, V_second, inflation_radius, t_halves.second);
    // Take the hull of the two halves of the vertices' trajectories
    for (int i = 0; i < vertices.rows(); i++) {
        for (int j = 0; j < vertices.cols(); j++) {
            vertices(i, j) = hull(V_first(i, j), V_second(i, j));
            assert(
                vertices(i, j).lower() - inflation_radius >= m_domainMin(j)
                && vertices(i, j).upper() + inflation_radius <= m_domainMax(j));
        }
    }

    return std::max(n_subs0, n_subs1) + 1;
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
        physics::cast<double, Interval>(poses_t1), /*body_ids,*/ vertices,
        inflation_radius);

    // Create a bounding box for all vertices
    std::vector<ipc::AABB> vertices_aabb;
    vertices_aabb.reserve(vertices.rows());
    for (long i = 0; i < vertices.rows(); i++) {
        Eigen::VectorX3d min(bodies.dim());
        Eigen::VectorX3d max(bodies.dim());
        for (int j = 0; j < vertices.cols(); j++) {
            assert(!empty(vertices(i, j)));
            min(j) = vertices(i, j).lower() - inflation_radius;
            max(j) = vertices(i, j).upper() + inflation_radius;
        }
        assert((min.array() >= m_domainMin.array()).all());
        assert((max.array() <= m_domainMax.array()).all());
        vertices_aabb.emplace_back(min, max);
    }

    // Add all vertices of the bodies
    for (long i = 0; i < vertices.rows(); i++) {
        this->addElement(vertices_aabb[i], i, this->m_vertexItems);
    }

    // Add all edge of the bodies
    for (long i = 0; i < bodies.m_edges.rows(); i++) {
        this->addElement(
            ipc::AABB(
                vertices_aabb[bodies.m_edges(i, 0)],
                vertices_aabb[bodies.m_edges(i, 1)]),
            i, this->m_edgeItems);
    }

    // Add all faces of the bodies
    for (long i = 0; i < bodies.m_faces.rows(); i++) {
        this->addElement(
            ipc::AABB(
                ipc::AABB(
                    vertices_aabb[bodies.m_faces(i, 0)],
                    vertices_aabb[bodies.m_faces(i, 1)]),
                vertices_aabb[bodies.m_faces(i, 2)]),
            i, this->m_faceItems);
    }
}

} // namespace ccd
