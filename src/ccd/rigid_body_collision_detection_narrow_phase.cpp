#include <ccd/rigid_body_collision_detection.hpp>

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
    const EdgeVertexCandidate& ev_candidate,
    EdgeVertexImpacts& ev_impacts)
{
    long bodyA_id = bodies.m_vertex_to_body_map(ev_candidate.vertex_index);
    const physics::RigidBody& bodyA = bodies.m_rbs[bodyA_id];
    long vertex_id =
        ev_candidate.vertex_index - bodies.m_body_vertex_id[bodyA_id];

    long bodyB_id = bodies.m_vertex_to_body_map(ev_candidate.vertex_index);
    const physics::RigidBody& bodyB = bodies.m_rbs[bodyB_id];
    long edge_id = ev_candidate.edge_index - bodies.m_body_edge_id[bodyB_id];

    double toi;
    bool are_colliding = compute_edge_vertex_time_of_impact(
        bodyA, poses[bodyA_id], displacements[bodyA_id], vertex_id, bodyB,
        poses[bodyB_id], displacements[bodyB_id], edge_id, toi);
    if (are_colliding) {
        ev_impacts.emplace_back(
            toi, ev_candidate.edge_index, /*alpha=*/-1,
            ev_candidate.vertex_index);
    }
}

void detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const EdgeEdgeCandidate& ee_candidate,
    EdgeEdgeImpacts& ee_impacts)
{
    throw NotImplementedError(
        "Edge-edge collisions not implemented for rigid bodies!");
}

void detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const FaceVertexCandidate& fv_candidate,
    FaceVertexImpacts& fv_impacts)
{
    throw NotImplementedError(
        "Face-vertex collisions not implemented for rigid bodies!");
}

} // namespace ccd
