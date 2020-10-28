#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
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

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_rigid_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::Candidates& candidates,
    ConcurrentImpacts& impacts);

/// @brief Determine if a single edge-vertext pair intersects.
bool edge_vertex_rigid_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& ev_candidate,
    double& toi,
    double& alpha,
    double earliest_toi = 1);

bool edge_edge_rigid_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& ee_candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    double earliest_toi = 1);

bool face_vertex_rigid_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& fv_candidate,
    double& toi,
    double& u,
    double& v,
    double earliest_toi = 1);

} // namespace ccd
