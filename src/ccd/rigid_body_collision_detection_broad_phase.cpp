#include "rigid_body_collision_detection.hpp"

#include <ccd/rigid_body_hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method)
{
    assert(bodies.num_bodies() == poses.size());
    assert(poses.size() == displacements.size());

    // Do the broad phase by detecting candidate impacts
    Candidates candidates;
    detect_collision_candidates(
        bodies, poses, displacements, collision_types, candidates, method);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_collisions_from_candidates(
        bodies, poses, displacements, candidates, impacts);
}

void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)
{
    ConcurrentImpacts impacts;
    detect_collisions(
        vertices, displacements, edges, Eigen::MatrixXi(), group_ids,
        CollisionType::EDGE_VERTEX, impacts, method);
    ev_impacts = impacts.ev_impacts;
}

void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)
{
    detect_edge_vertex_collisions(
        vertices, displacements, edges, Eigen::VectorXi(), ev_impacts, method);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    Candidates& candidates,
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
            bodies.m_vertex_to_body_map, collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_hash_grid(
            bodies, poses, displacements, collision_types, candidates,
            inflation_radius);
        spdlog::debug(
            "hash_grid_ev_candidates.size()={:d} "
            "hash_grid_ee_candidates.size()={:d} "
            "hash_grid_fv_candidates.size()={:d}",
            candidates.ev_candidates.size(), candidates.ee_candidates.size(),
            candidates.fv_candidates.size());
        break;
    }

    PROFILE_END(BROAD_PHASE);
    PROFILE_END();
}

void detect_edge_vertex_collision_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    Candidates candidates;
    detect_collision_candidates(
        vertices, displacements, edges, Eigen::MatrixXi(), group_ids,
        CollisionType::EDGE_VERTEX, candidates, method, inflation_radius);
    ev_candidates = candidates.ev_candidates;
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    Candidates& candidates,
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
                    bodies.m_edges, group_ids, candidates.ev_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::EDGE_EDGE) {
                hashgrid.getEdgeEdgePairs(
                    bodies.m_edges, group_ids, candidates.ee_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::FACE_VERTEX) {
                hashgrid.getFaceVertexPairs(
                    bodies.m_faces, group_ids, candidates.fv_candidates);
            }
        });
}

} // namespace ccd
