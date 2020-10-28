#pragma once

#include <Eigen/Core>

#include <nlohmann/json.hpp>

#include <ipc/spatial_hash/collision_candidate.hpp>

#include <ccd/impact.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

/// @brief Possible methods for detecting all edge vertex collisions.
enum DetectionMethod {
    BRUTE_FORCE, ///< @brief Use brute-force to detect all collisions.
    HASH_GRID ///< @brief Use a spatial data structure to detect all collisions.
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    DetectionMethod,
    { { HASH_GRID, "hash_grid" }, { BRUTE_FORCE, "brute_force" } });

/// @brief Possible trajectories of vertices in a rigid body.
enum TrajectoryType {
    /// @brief Linearization of the rotation component of rigid body
    /// trajectories.
    LINEAR,
    /// @brief Piceiwise linearization of the the rotation component of rigid
    /// body trajectories.
    PIECEWISE_LINEAR,
    /// @brief Fully nonlinear rigid trajectories.
    RIGID
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    TrajectoryType,
    { { LINEAR, "linear" },
      { PIECEWISE_LINEAR, "piecewise_linear" },
      { RIGID, "rigid" } });

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
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    TrajectoryType trajectory = TrajectoryType::RIGID);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    TrajectoryType trajectory = TrajectoryType::RIGID,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::Candidates& candidates,
    ConcurrentImpacts& impacts,
    TrajectoryType trajectory = TrajectoryType::RIGID);

/// @brief Determine if a single edge-vertext pair intersects.
bool edge_vertex_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& ev_candidate,
    double& toi,
    double& alpha,
    TrajectoryType trajectory = TrajectoryType::RIGID,
    double earliest_toi = 1);

bool edge_edge_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& ee_candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    TrajectoryType trajectory = TrajectoryType::RIGID,
    double earliest_toi = 1);

bool face_vertex_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& fv_candidate,
    double& toi,
    double& u,
    double& v,
    TrajectoryType trajectory = TrajectoryType::RIGID,
    double earliest_toi = 1);

} // namespace ccd
