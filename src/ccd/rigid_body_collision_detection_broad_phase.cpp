#include <ccd/rigid_body_collision_detection.hpp>

#include <ccd/rigid_body_hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const int collision_types,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts,
    DetectionMethod method)
{
    assert(bodies.num_bodies() == poses.size());
    assert(poses.size() == displacements.size());

    // Do the broad phase by detecting candidate impacts
    EdgeVertexCandidates ev_candidates;
    EdgeEdgeCandidates ee_candidates;
    FaceVertexCandidates fv_candidates;
    detect_collision_candidates(
        bodies, poses, displacements, collision_types, ev_candidates,
        ee_candidates, fv_candidates, method);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_collisions_from_candidates(
        bodies, poses, displacements, ev_candidates, ee_candidates,
        fv_candidates, ev_impacts, ee_impacts, fv_impacts);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    assert(
        method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);

    PROFILE_START();
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            bodies.world_vertices(poses), bodies.m_edges, bodies.m_faces,
            bodies.m_vertex_to_body_map, collision_types, ev_candidates,
            ee_candidates, fv_candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_hash_grid(
            bodies, poses, displacements, collision_types, ev_candidates,
            ee_candidates, fv_candidates, inflation_radius);
        spdlog::debug(
            "hash_grid_ev_candidates.size()={:d} "
            "hash_grid_ee_candidates.size()={:d} "
            "hash_grid_fv_candidates.size()={:d}",
            ev_candidates.size(), ee_candidates.size(), fv_candidates.size());
        break;
    }

    PROFILE_END(BROAD_PHASE);
    PROFILE_END();
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<physics::Pose<double>>& poses,
    const std::vector<physics::Pose<double>>& displacements,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    const double inflation_radius)
{
    RigidBodyHashGrid hashgrid;
    hashgrid.resize(bodies, poses, displacements, inflation_radius);
    hashgrid.addBodies(bodies, poses, displacements, inflation_radius);

    const Eigen::VectorXi& group_ids = bodies.m_vertex_to_body_map;
    tbb::parallel_invoke(
        [&] {
            if (collision_types & CollisionType::EDGE_VERTEX) {
                hashgrid.getVertexEdgePairs(
                    bodies.m_edges, group_ids, ev_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::EDGE_EDGE) {
                hashgrid.getEdgeEdgePairs(
                    bodies.m_edges, group_ids, ee_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::FACE_VERTEX) {
                hashgrid.getFaceVertexPairs(
                    bodies.m_faces, group_ids, fv_candidates);
            }
        });
}

} // namespace ccd
