#include "rigid_body_collision_detection.hpp"

#include <cmath>

#include <ccd/rigid_body_time_of_impact.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const Candidates& candidates,
    ConcurrentImpacts& impacts)
{
#ifdef LINEARIZED_TRAJECTORY_CCD
    Eigen::MatrixXd V_t0 = bodies.world_vertices(poses);
    Eigen::MatrixXd V_t1 = bodies.world_vertices(poses + displacements);
    detect_collisions_from_candidates(
        V_t0, V_t1 - V_t0, bodies.m_edges, bodies.m_faces, candidates, impacts);
#else
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    auto detect_ev_collision = [&](const EdgeVertexCandidate& ev_candidate) {
        double toi, alpha;
        bool are_colliding = detect_edge_vertex_collisions_narrow_phase(
            bodies, poses, displacements, ev_candidate, toi, alpha);
        if (are_colliding) {
            impacts.ev_impacts.emplace_back(
                toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index);
        }
    };

    auto detect_ee_collision = [&](const EdgeEdgeCandidate& ee_candidate) {
        double toi, edge0_alpha, edge1_alpha;
        bool are_colliding = detect_edge_edge_collisions_narrow_phase(
            bodies, poses, displacements, ee_candidate, toi, edge0_alpha,
            edge1_alpha);
        if (are_colliding) {
            impacts.ee_impacts.emplace_back(
                toi, ee_candidate.edge0_index, edge0_alpha,
                ee_candidate.edge1_index, edge1_alpha);
        }
    };

    auto detect_fv_collision = [&](const FaceVertexCandidate& fv_candidate) {
        double toi, u, v;
        bool are_colliding = detect_face_vertex_collisions_narrow_phase(
            bodies, poses, displacements, fv_candidate, toi, u, v);
        if (are_colliding) {
            impacts.fv_impacts.emplace_back(
                toi, fv_candidate.face_index, u, v, fv_candidate.vertex_index);
        }
    };

    tbb::parallel_invoke(
        [&] {
            impacts.ev_impacts.clear();
            tbb::parallel_for_each(
                candidates.ev_candidates, detect_ev_collision);
        },

        [&] {
            impacts.ee_impacts.clear();
            tbb::parallel_for_each(
                candidates.ee_candidates, detect_ee_collision);
        },

        [&] {
            impacts.fv_impacts.clear();
            tbb::parallel_for_each(
                candidates.fv_candidates, detect_fv_collision);
        });

    PROFILE_END(NARROW_PHASE);
    PROFILE_END();
#endif
}

// Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeVertexCandidate& candidate,
    double& toi,
    double& alpha)
{
#ifdef LINEARIZED_TRAJECTORY_CCD
    Eigen::MatrixXd V_t0 = bodies.world_vertices(poses);
    Eigen::MatrixXd V_t1 = bodies.world_vertices(poses + displacements);
    Eigen::MatrixXd U = V_t1 - V_t0;
    return detect_edge_vertex_collisions_narrow_phase(
        V_t0.row(bodies.m_edges(candidate.edge_index, 0)),
        V_t0.row(bodies.m_edges(candidate.edge_index, 1)),
        V_t0.row(candidate.vertex_index),
        U.row(bodies.m_edges(candidate.edge_index, 0)),
        U.row(bodies.m_edges(candidate.edge_index, 1)),
        U.row(candidate.vertex_index), toi, alpha);
#else
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    bool are_colliding = compute_edge_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edge_id, toi);
    if (are_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA =
            poses[bodyA_id] + displacements[bodyA_id] * toi;
        physics::Pose<double> poseB =
            poses[bodyB_id] + displacements[bodyB_id] * toi;

        // Get the world vertex of the point at time t
        Eigen::VectorX3d v0 =
            bodies.m_rbs[bodyA_id].world_vertex(poseA, vertex_id);
        // Get the world vertex of the edge at time t
        Eigen::VectorX3d v1 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 0));
        Eigen::VectorX3d v2 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 1));

        // Project the point onto the edge by computing its scalar projection
        Eigen::VectorX3d edge_vec = v2 - v1;
        alpha = (v0 - v1).dot(edge_vec) / edge_vec.squaredNorm();
    }
    return are_colliding;
#endif
}

bool detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeEdgeCandidate& candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha)
{
#ifdef LINEARIZED_TRAJECTORY_CCD
    throw NotImplementedError("Conversion from RB to linearized trajectories");
#else
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);

    bool are_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        edgeA_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edgeB_id, toi);
    if (are_colliding) {
        edge0_alpha = -1; // TODO: Compute this correctly
        edge0_alpha = -1; // TODO: Compute this correctly
    }
    return are_colliding;
#endif
}

bool detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const FaceVertexCandidate& candidate,
    double& toi,
    double& u,
    double& v)
{
#ifdef LINEARIZED_TRAJECTORY_CCD
    throw NotImplementedError("Conversion from RB to linearized trajectories");
#else
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);

    bool are_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], face_id, toi);
    if (are_colliding) {
        u = -1; // TODO: Compute this correctly
        v = -1; // TODO: Compute this correctly
    }
    return are_colliding;
#endif
}

} // namespace ccd