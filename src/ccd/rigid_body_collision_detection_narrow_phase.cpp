#include "rigid_body_collision_detection.hpp"

#include <ccd/rigid_body_time_of_impact.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const EdgeVertexCandidates& ev_candidates,
    const EdgeEdgeCandidates& ee_candidates,
    const FaceVertexCandidates& fv_candidates,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts)
{
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    tbb::parallel_invoke(
        [&] {
            ev_impacts.clear();
            tbb::parallel_for_each(
                ev_candidates, [&](const EdgeVertexCandidate& ev_candidate) {
                    detect_edge_vertex_collisions_narrow_phase(
                        bodies, poses, displacements, ev_candidate, ev_impacts);
                });
        },

        [&] {
            ee_impacts.clear();
            tbb::parallel_for_each(
                ee_candidates, [&](const EdgeEdgeCandidate& ee_candidate) {
                    detect_edge_edge_collisions_narrow_phase(
                        bodies, poses, displacements, ee_candidate, ee_impacts);
                });
        },

        [&] {
            fv_impacts.clear();
            tbb::parallel_for_each(
                fv_candidates, [&](const FaceVertexCandidate& fv_candidate) {
                    detect_face_vertex_collisions_narrow_phase(
                        bodies, poses, displacements, fv_candidate, fv_impacts);
                });
        });

    PROFILE_END(NARROW_PHASE);
    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
void detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const EdgeVertexCandidate& candidate,
    EdgeVertexImpacts& impacts)
{
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    double toi;
    bool are_colliding = compute_edge_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edge_id, toi);
    if (are_colliding) {
        impacts.emplace_back(
            toi, candidate.edge_index, /*alpha=*/-1, candidate.vertex_index);
    }
}

void detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const EdgeEdgeCandidate& candidate,
    EdgeEdgeImpacts& impacts)
{
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);

    double toi;
    bool are_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        edgeA_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edgeB_id, toi);
    if (are_colliding) {
        impacts.emplace_back(
            toi, candidate.edge0_index, /*impacted_alpha=*/-1,
            candidate.edge1_index, /*impacting_alpha=*/-1);
    }
}

void detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const FaceVertexCandidate& candidate,
    FaceVertexImpacts& impacts)
{
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);

    double toi;
    bool are_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], face_id, toi);
    if (are_colliding) {
        impacts.emplace_back(
            toi, candidate.face_index, /*u=*/-1, /*v=*/-1,
            candidate.vertex_index);
    }
}

} // namespace ccd
