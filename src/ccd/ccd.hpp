#pragma once

#include <Eigen/Core>

#include <nlohmann/json.hpp>

#include <ipc/broad_phase/collision_candidate.hpp>

#include <ccd/impact.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

/// @brief Possible methods for detecting all edge vertex collisions.
enum DetectionMethod {
    BRUTE_FORCE, ///< @brief Use brute-force to detect all collisions
    HASH_GRID, ///< @brief Use a spatial data structure to detect all collisions
    BVH,       ///< @brief Use a BVH to detect all collisions
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    DetectionMethod,
    { { HASH_GRID, "hash_grid" },
      { BRUTE_FORCE, "brute_force" },
      { BVH, "bvh" } });

/// @brief Possible trajectories of vertices in a rigid body.
enum TrajectoryType {
    /// @brief Linearization of the rotation component of rigid body
    /// trajectories.
    LINEAR,
    /// @brief Piceiwise linearization of the the rotation component of rigid
    /// body trajectories. CCD is computed conservativly using minimum
    /// separation CCD to bound the difference between RIGID.
    PIECEWISE_LINEAR,
    /// @brief Fully nonlinear rigid trajectories.
    RIGID,
    /// @brief Same trajectory as RIGID, but the time of impact is computed
    /// using Redon et al. [2002].
    REDON
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    TrajectoryType,
    { { LINEAR, "linear" },
      { PIECEWISE_LINEAR, "piecewise_linear" },
      { RIGID, "rigid" },
      { REDON, "redon" } });

namespace CollisionType {
    static const int EDGE_VERTEX = 1;
    static const int EDGE_EDGE = 2;
    static const int FACE_VERTEX = 4;
} // namespace CollisionType

///////////////////////////////////////////////////////////////////////////////
// CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Find all collisions in one time step.
void detect_collisions(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Impacts& impacts,
    DetectionMethod method,
    TrajectoryType trajectory);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    TrajectoryType trajectory,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const Candidates& candidates,
    Impacts& impacts,
    TrajectoryType trajectory);

/// @brief Determine if a single edge-vertext pair intersects.
bool edge_vertex_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& ev_candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi = 1,
    double minimum_separation_distance = 0);

bool edge_edge_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& ee_candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi = 1,
    double minimum_separation_distance = 0);

bool face_vertex_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& fv_candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi = 1,
    double minimum_separation_distance = 0);

double edge_vertex_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& candidate,
    double toi,
    TrajectoryType trajectory);

void edge_edge_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& candidate,
    double toi,
    double& alpha,
    double& beta,
    TrajectoryType trajectory);

void face_vertex_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& candidate,
    double toi,
    double& u,
    double& v,
    TrajectoryType trajectory);

} // namespace ipc::rigid
