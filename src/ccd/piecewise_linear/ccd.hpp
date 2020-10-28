#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Determine if a single edge-vertext pair intersects.
bool edge_vertex_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& ev_candidate,
    double& toi,
    double& alpha,
    double earliest_toi = 1);

bool edge_edge_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& ee_candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    double earliest_toi = 1);

bool face_vertex_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& fv_candidate,
    double& toi,
    double& u,
    double& v,
    double earliest_toi = 1);

} // namespace ccd
