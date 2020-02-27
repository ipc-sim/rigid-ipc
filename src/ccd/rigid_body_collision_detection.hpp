#pragma once

#include <vector>

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <physics/rigid_body_assembler.hpp>

// #define LINEARIZED_TRAJECTORY_CCD

namespace ccd {

/// @brief Find all collisions in one time step.
void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID);

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
    ConcurrentImpacts& impacts);

/// @brief Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeVertexCandidate& ev_candidate,
    double& toi,
    double& alpha);

bool detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeEdgeCandidate& ee_candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha);

bool detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const FaceVertexCandidate& fv_candidate,
    double& toi,
    double& u,
    double& v);

} // namespace ccd
