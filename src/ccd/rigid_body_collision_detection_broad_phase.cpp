#include "rigid_body_collision_detection.hpp"

#include <ccd/rigid_body_hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method,
    TrajectoryType trajectory)
{
    assert(bodies.num_bodies() == poses_t0.size());
    assert(poses_t0.size() == poses_t1.size());

    // Do the broad phase by detecting candidate impacts
    Candidates candidates;
    detect_collision_candidates(
        bodies, poses_t0, poses_t1, collision_types, candidates, method,
        trajectory);

    // Do the narrow phase by detecting actual impacts from the candidate
    // set
    detect_collisions_from_candidates(
        bodies, poses_t0, poses_t1, candidates, impacts, trajectory);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    TrajectoryType trajectory,
    const double inflation_radius)
{
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    // TODO: Add a quick capsule test for intersecting rigid bodies
    // def capsule(body, pose, displacement):
    //     return {"radius": body.r_max, "start": pose.position, "end":
    //     pose.position + displacement.position}
    //
    // def distance(capsule1, capsule2):
    //     return segment_segment_distance(capsule1.start, capsule1.end,
    //     capsule2.start, capsule2.end) - capsule1.radius - capsule2.radius
    //
    // capsules = [capsule(body, pose, displacement) for body in bodies]
    //
    // intersecting_capsules = {}
    // # O(n^2) segment_segment_distance calls
    // for i in range(len(capsules)):
    //     for j in range(i, len(capsules)):
    //         if distance(capsules[i], capsules[j]) <= tol:
    //             intersecting_capsules.add(i)
    //             intersecting_capsules.add(j)
    //
    // detect_collisions(bodies[intersecting_capsules], ...)

    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses_t0);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses_t1);
        detect_collision_candidates(
            V_t0, V_t1, bodies.m_edges, bodies.m_faces, bodies.group_ids(),
            collision_types, candidates, method, inflation_radius);
        return;
    }

    assert(trajectory == TrajectoryType::SCREWING);
    assert(
        method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);

    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            bodies.world_vertices(poses_t0), bodies.m_edges, bodies.m_faces,
            bodies.group_ids(), collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_hash_grid(
            bodies, poses_t0, poses_t1, collision_types, candidates,
            inflation_radius);
        break;
    }

    PROFILE_END(BROAD_PHASE);
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius)
{
    RigidBodyHashGrid hashgrid;
    hashgrid.resize(bodies, poses_t0, poses_t1, inflation_radius);
    hashgrid.addBodies(bodies, poses_t0, poses_t1, inflation_radius);

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
