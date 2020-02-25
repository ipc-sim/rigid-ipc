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
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    bool are_colliding = compute_edge_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edge_id, toi);
    if (are_colliding) {
        alpha = -1; // TODO: Compute this correctly
    }
    return are_colliding;
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
}

} // namespace ccd
