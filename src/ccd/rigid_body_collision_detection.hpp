#pragma once

#include <vector>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <ccd/collision_detection.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

enum TrajectoryType {
    LINEARIZED, ///< Linearize the rotation component of rigid body trajectories
    SCREWING    ///< Use the fully nonlinear screwing trajectories
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    TrajectoryType, { { LINEARIZED, "linearized" }, { SCREWING, "screwing" } });

/// @brief Find all collisions in one time step.
void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    TrajectoryType trajectory = TrajectoryType::SCREWING);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    TrajectoryType trajectory = TrajectoryType::SCREWING,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const Candidates& candidates,
    ConcurrentImpacts& impacts,
    TrajectoryType trajectory = TrajectoryType::SCREWING);

/// @brief Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeVertexCandidate& ev_candidate,
    double& toi,
    double& alpha,
    TrajectoryType trajectory = TrajectoryType::SCREWING);

bool detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeEdgeCandidate& ee_candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    TrajectoryType trajectory = TrajectoryType::SCREWING);

bool detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const FaceVertexCandidate& fv_candidate,
    double& toi,
    double& u,
    double& v,
    TrajectoryType trajectory = TrajectoryType::SCREWING);

} // namespace ccd
