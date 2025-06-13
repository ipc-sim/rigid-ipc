#include "rigid_body_assembler.hpp"

#include <Eigen/Geometry>
#include <finitediff.hpp>
#include <igl/PI.h>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <tbb/parallel_sort.h>

#include <logger.hpp>
#include <physics/mass.hpp>
#include <profiler.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/flatten.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

void RigidBodyAssembler::init(const std::vector<RigidBody>& rigid_bodies)
{
    m_rbs = rigid_bodies;

    size_t num_bodies = rigid_bodies.size();
    m_body_vertex_id.resize(num_bodies + 1);
    m_body_face_id.resize(num_bodies + 1);
    m_body_edge_id.resize(num_bodies + 1);

    // Store the starting position of each RB vertex in the global vertices
    m_body_vertex_id[0] = m_body_face_id[0] = m_body_edge_id[0] = 0;
    for (size_t i = 0; i < num_bodies; ++i) {
        auto& rb = rigid_bodies[i];
        m_body_vertex_id[i + 1] = m_body_vertex_id[i] + rb.num_vertices();
        m_body_face_id[i + 1] = m_body_face_id[i] + rb.faces.rows();
        m_body_edge_id[i + 1] = m_body_edge_id[i] + rb.edges.rows();
    }

    // global edges and faces
    m_edges.resize(m_body_edge_id.back(), 2);
    m_faces.resize(m_body_face_id.back(), 3);
    m_faces_to_edges.resize(m_body_face_id.back(), 3);
    m_codim_edges_to_edges.clear();
    m_codim_vertices_to_vertices.clear();
    for (size_t i = 0; i < num_bodies; ++i) {
        auto& rb = rigid_bodies[i];
        if (rb.edges.size() != 0) {
            m_edges.block(m_body_edge_id[i], 0, rb.edges.rows(), 2) =
                rb.edges.array() + m_body_vertex_id[i];
        }
        if (rb.faces.size() != 0) {
            m_faces.block(m_body_face_id[i], 0, rb.faces.rows(), 3) =
                rb.faces.array() + m_body_vertex_id[i];
            m_faces_to_edges.block(m_body_face_id[i], 0, rb.faces.rows(), 3) =
                rb.mesh_selector.face_to_edges().array() + m_body_edge_id[i];
        }
        for (const auto& ei : rb.mesh_selector.codim_edges_to_edges()) {
            m_codim_edges_to_edges.push_back(ei + m_body_edge_id[i]);
        }
        for (const auto& vi : rb.mesh_selector.codim_vertices_to_vertices()) {
            m_codim_vertices_to_vertices.push_back(vi + m_body_vertex_id[i]);
        }
    }
    // vertex to body map
    m_vertex_to_body_map.resize(num_vertices());
    for (size_t i = 0; i < num_bodies; ++i) {
        auto& rb = rigid_bodies[i];
        m_vertex_to_body_map.segment(m_body_vertex_id[i], rb.num_vertices())
            .setConstant(int(i));
    }
    // vertex to group id map
    m_vertex_group_ids.resize(num_vertices());
    for (size_t i = 0; i < num_bodies; ++i) {
        auto& rb = rigid_bodies[i];
        m_vertex_group_ids.segment(m_body_vertex_id[i], rb.num_vertices())
            .setConstant(rb.group_id);
    }

    // rigid body mass-matrix
    int rb_ndof = num_bodies ? rigid_bodies[0].ndof() : 0;
    m_rb_mass_matrix.resize(num_bodies * rb_ndof);
    for (int i = 0; i < int(num_bodies); ++i) {
        m_rb_mass_matrix.diagonal().segment(i * rb_ndof, rb_ndof) =
            rigid_bodies[i].mass_matrix.diagonal();
    }

    // rigid_body dof_fixed flag
    is_rb_dof_fixed.resize(num_bodies * rb_ndof);
    for (int i = 0; i < int(num_bodies); ++i) {
        auto& rb = rigid_bodies[size_t(i)];
        is_rb_dof_fixed.segment(rb_ndof * i, rb_ndof) = rb.is_dof_fixed;
    }

    // rigid_body vertex dof_fixed flag
    is_dof_fixed.resize(num_vertices(), rb_ndof);
    for (size_t i = 0; i < num_bodies; ++i) {
        auto& rb = rigid_bodies[i];
        is_dof_fixed.block(m_body_vertex_id[i], 0, rb.num_vertices(), rb_ndof) =
            rb.is_dof_fixed.transpose().replicate(rb.num_vertices(), 1);
    }

    average_edge_length = 0;
    for (const auto& body : rigid_bodies) {
        average_edge_length += body.edges.rows() * body.average_edge_length;
    }
    average_edge_length /= m_edges.rows();
    assert(std::isfinite(average_edge_length));

    average_mass = 0;
    int num_free_dof = 0;
    for (int i = 0; i < is_rb_dof_fixed.size(); i++) {
        if (!is_rb_dof_fixed[i]) {
            average_mass += m_rb_mass_matrix.diagonal()[i];
            num_free_dof++;
        }
    }
    if (num_free_dof) {
        average_mass /= num_free_dof;
    }
}

size_t RigidBodyAssembler::count_kinematic_bodies() const
{
    size_t n = 0;
    for (const auto& body : m_rbs) {
        if (body.type == RigidBodyType::KINEMATIC) {
            n++;
        }
    }
    return n;
}

void RigidBodyAssembler::global_to_local_vertex(
    const long global_vertex_id,
    long& rigid_body_id,
    long& local_vertex_id) const
{
    rigid_body_id = vertex_id_to_body_id(global_vertex_id);
    local_vertex_id =
        global_vertex_id - m_body_vertex_id[size_t(rigid_body_id)];
}

void RigidBodyAssembler::global_to_local_edge(
    const long global_edge_id, long& rigid_body_id, long& local_edge_id) const
{
    rigid_body_id = edge_id_to_body_id(global_edge_id);
    local_edge_id = global_edge_id - m_body_edge_id[size_t(rigid_body_id)];
}

void RigidBodyAssembler::global_to_local_face(
    const long global_face_id, long& rigid_body_id, long& local_face_id) const
{
    rigid_body_id = face_id_to_body_id(global_face_id);
    local_face_id = global_face_id - m_body_face_id[size_t(rigid_body_id)];
}

PosesD RigidBodyAssembler::rb_poses(const bool previous) const
{
    PosesD poses;
    poses.resize(num_bodies());
    for (size_t i = 0; i < num_bodies(); i++) {
        poses[i] = previous ? m_rbs[i].pose_prev : m_rbs[i].pose;
    }
    return poses;
}

void RigidBodyAssembler::set_rb_poses(const PosesD& poses)
{
    assert(num_bodies() == poses.size());
    for (size_t i = 0; i < num_bodies(); i++) {
        m_rbs[i].pose = poses[i];
        // if (m_rbs[i].pose.rotation.size() > 1) {
        //     // Mod 2π on the angle
        //     double angle = m_rbs[i].pose.rotation.norm();
        //     if (abs(angle) >= 2 * igl::PI) {
        //         m_rbs[i].pose.rotation /= angle;
        //         m_rbs[i].pose.rotation *= fmod(angle, 2 * igl::PI);
        //     }
        // }
    }
}

Eigen::MatrixXd
RigidBodyAssembler::world_vertices(const RigidBody::Step step) const
{
    Eigen::MatrixXd V(num_vertices(), dim());
    for (size_t i = 0; i < num_bodies(); ++i) {
        auto& rb = m_rbs[i];
        V.block(m_body_vertex_id[i], 0, rb.num_vertices(), dim()) =
            rb.world_vertices(step);
    }
    return V;
}

Eigen::MatrixXd RigidBodyAssembler::world_velocities() const
{
    Eigen::MatrixXd V(num_vertices(), dim());
    for (size_t i = 0; i < num_bodies(); ++i) {
        auto& rb = m_rbs[i];
        V.block(m_body_vertex_id[i], 0, rb.num_vertices(), dim()) =
            rb.world_velocities();
    }
    return V;
}

Eigen::MatrixXd RigidBodyAssembler::world_vertices_diff(
    const PosesD& poses,
    Eigen::MatrixXd& jac,
    Eigen::MatrixXd& hess,
    bool compute_jac,
    bool compute_hess) const
{
    assert(num_bodies() == poses.size());

    if (!compute_jac && !compute_hess) {
        return world_vertices(poses);
    }

    // We will only use auto diff to compute the derivatives of the rotation
    // matrix.
    typedef AutodiffType<Eigen::Dynamic, /*maxN=*/3> Diff;

    PROFILE_POINT("RigidBodyAssembler::world_vertices_diff");
    PROFILE_START();

    // V: Rᵐ ↦ Rⁿ (vertices flattened rowwise)
    const int n = num_vertices() * dim();
    const int m = PoseD::dim_to_ndof(dim());

    NAMED_PROFILE_POINT(
        "RigidBodyAssembler::world_vertices_diff:allocation", ALLOCATION);
    PROFILE_START(ALLOCATION);
    if (compute_jac) {
        // ∇ V(x): Rᵐ ↦ Rⁿˣᵐ
        jac.resize(n, m);
    }
    if (compute_hess) {
        // ∇²V(x): Rᵐ ↦ Rⁿˣᵐˣᵐ
        hess.resize(n * m, m); // Store as (n * m) × m matrix
    }
    Eigen::MatrixXd V(num_vertices(), dim());
    PROFILE_END(ALLOCATION);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), num_bodies()),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t rb_i = range.begin(); rb_i != range.end(); ++rb_i) {
                const RigidBody& rb = m_rbs[rb_i];

                // Index of rigid bodies first vertex in the global vertices
                long rb_v0_i = m_body_vertex_id[rb_i];

                if (compute_hess) {
                    rb.world_vertices_diff<Diff::DDouble2>(
                        poses[rb_i], rb_v0_i, V, jac, hess);
                } else {
                    assert(compute_jac);
                    rb.world_vertices_diff<Diff::DDouble1>(
                        poses[rb_i], rb_v0_i, V, jac, hess);
                }
            }
        });

    assert((V - world_vertices(poses)).norm() < 1e-12);

    PROFILE_END();
    return V;
}

std::vector<std::pair<int, int>> RigidBodyAssembler::close_bodies(
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const double inflation_radius) const
{
    // if (num_bodies() < 10) {
    //     return close_bodies_brute_force(
    //         poses_t0, poses_t1, inflation_radius);
    // }
    //
    // return close_bodies_hash_grid(poses_t0, poses_t1, inflation_radius);
    return close_bodies_bvh(poses_t0, poses_t1, inflation_radius);
}

std::vector<std::pair<int, int>> RigidBodyAssembler::close_bodies_brute_force(
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const double inflation_radius) const
{
    std::vector<std::pair<int, int>> close_body_pairs;
    for (int i = 0; i < num_bodies(); i++) {
        double ri = m_rbs[i].r_max;
        for (int j = i + 1; j < num_bodies(); j++) {
            if (m_rbs[i].group_id == m_rbs[j].group_id) {
                continue;
            }

            double rj = m_rbs[j].r_max;
            double distance = sqrt(edge_edge_distance(
                poses_t0[i].position, poses_t1[i].position,
                poses_t0[j].position, poses_t1[j].position));
            if (distance <= ri + rj + inflation_radius) {
                close_body_pairs.emplace_back(i, j);
            }
        }
    }
    return close_body_pairs;
}

std::vector<std::pair<int, int>> RigidBodyAssembler::close_bodies_bvh(
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const double inflation_radius) const
{
    std::vector<std::array<Eigen::Vector3d, 2>> body_bounding_boxes;
    body_bounding_boxes.reserve(num_bodies());

    NAMED_PROFILE_POINT("RigidBodyAssembler::close_bodies_bvh:build", BUILD);
    PROFILE_START(BUILD);

    for (int i = 0; i < num_bodies(); i++) {
        VectorMax3d min, max;
        m_rbs[i].compute_bounding_box(poses_t0[i], poses_t1[i], min, max);
        min.array() -= inflation_radius;
        max.array() += inflation_radius;
        Eigen::Vector3d min3D = Eigen::Vector3d::Zero(),
                        max3D = Eigen::Vector3d::Zero();
        min3D.head(dim()) = min;
        max3D.head(dim()) = max;
        body_bounding_boxes.push_back({ { min3D, max3D } });
    }

    BVH::BVH bvh;
    bvh.init(body_bounding_boxes);

    PROFILE_END(BUILD);

    NAMED_PROFILE_POINT("RigidBodyAssembler::close_bodies_bvh:query", QUERY);
    PROFILE_START(QUERY);

    std::vector<std::pair<int, int>> close_body_pairs;
    for (int i = 0; i < num_bodies(); i++) {
        std::vector<unsigned int> intersecting_body_ids;
        bvh.intersect_box(
            body_bounding_boxes[i][0], body_bounding_boxes[i][1],
            intersecting_body_ids);
        for (const auto& j : intersecting_body_ids) {
            if (i < j && m_rbs[i].group_id != m_rbs[j].group_id) {
                close_body_pairs.emplace_back(i, j);
            }
        }
    }

    PROFILE_END(QUERY);
    PROFILE_MESSAGE(
        QUERY, "num_pairs", fmt::format("{:d}", close_body_pairs.size()));

    return close_body_pairs;
}

std::vector<std::pair<int, int>> RigidBodyAssembler::close_bodies_hash_grid(
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const double inflation_radius) const
{
    HashGrid hashgrid;

    Eigen::MatrixXd V(2 * num_bodies(), dim());
    Eigen::MatrixXi E(num_bodies(), dim());
    Eigen::VectorXi group_ids(2 * num_bodies());
    double max_radius = 0;
    for (int i = 0; i < num_bodies(); i++) {
        V.row(2 * i + 0) = poses_t0[i].position;
        V.row(2 * i + 1) = poses_t1[i].position;
        E(i, 0) = 2 * i + 0;
        E(i, 1) = 2 * i + 1;
        max_radius = std::max(m_rbs[i].r_max, max_radius);
        group_ids[2 * i + 1] = group_ids[2 * i] = m_rbs[i].group_id;
    }

    // Resize the grid
    hashgrid.resize(V, V, E, max_radius + inflation_radius);

    // Add each body
    for (int i = 0; i < num_bodies(); i++) {
        hashgrid.addEdge(
            V.row(E(i, 0)), V.row(E(i, 1)), V.row(E(i, 0)), V.row(E(i, 1)), i,
            m_rbs[i].r_max);
    }

    auto can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids[vi] != group_ids[vj];
    };

    std::vector<EdgeEdgeCandidate> body_candidates;
    hashgrid.getEdgeEdgePairs(E, body_candidates, can_collide);

    std::vector<std::pair<int, int>> close_body_pairs;
    close_body_pairs.reserve(body_candidates.size());
    for (const auto& body_candidate : body_candidates) {
        close_body_pairs.emplace_back(
            body_candidate.edge0_index, body_candidate.edge1_index);
    }

    return close_body_pairs;
}

} // namespace ipc::rigid
