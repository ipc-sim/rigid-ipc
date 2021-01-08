#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Discrete Collision Detection
// NOTE: Yes, I know this is inside the CCD directory.
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Continous Collision Detection
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

} // namespace ccd
