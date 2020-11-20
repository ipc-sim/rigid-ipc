#include "broad_phase.hpp"

#include <tbb/parallel_invoke.h>

#include <ccd/linear/broad_phase.hpp>
#include <ccd/rigid/rigid_body_hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    assert(method == BRUTE_FORCE || method == HASH_GRID);
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            bodies.world_vertices(poses_t0), bodies.m_edges, bodies.m_faces,
            bodies.group_ids(), collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_rigid_hash_grid(
            bodies, poses_t0, poses_t1, collision_types, candidates,
            inflation_radius);
        break;
    }

    PROFILE_END(BROAD_PHASE);
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    std::vector<int> body_ids =
        bodies.close_bodies(poses_t0, poses_t1, inflation_radius);

    if (body_ids.size() <= 1) {
        return;
    }

    RigidBodyHashGrid hashgrid;
    hashgrid.resize(bodies, poses_t0, poses_t1, body_ids, inflation_radius);
    hashgrid.addBodies(bodies, poses_t0, poses_t1, body_ids, inflation_radius);

    const Eigen::VectorXi& group_ids = bodies.group_ids();
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
