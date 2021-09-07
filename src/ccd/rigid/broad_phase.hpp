#pragma once

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/broad_phase/hash_grid.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <interval/interval.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Discrete Collision Detection
// NOTE: Yes, I know this is inside the CCD directory.
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const RigidBodyAssembler& bodies,
    const PosesD& poses,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const RigidBodyAssembler& bodies,
    const PosesD& poses,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Continous Collision Detection
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Intersection Detection
///////////////////////////////////////////////////////////////////////////////

/// Use a BVH to create a set of all candidate intersections.
void detect_intersection_candidates_rigid_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses,
    std::vector<EdgeFaceCandidate>& ef_candidates);

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

typedef tbb::enumerable_thread_specific<Candidates>
    ThreadSpecificCandidates;

void merge_local_candidates(
    const ThreadSpecificCandidates& storages, Candidates& candidates);

} // namespace ipc::rigid
