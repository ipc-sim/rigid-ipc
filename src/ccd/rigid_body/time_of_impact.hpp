/// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <physics/rigid_body.hpp>

/**
 * @namespace ccd:
 * @brief
 */
namespace ccd {

/// Find time-of-impact between two rigid bodies
bool compute_edge_vertex_time_of_impact(
    const physics::RigidBody& bodyA,
    const physics::Pose<double>& displacementA, // Displacement of bodyA
    const size_t& vertex_id,                    // In bodyA
    const physics::RigidBody& bodyB,
    const physics::Pose<double>& displacementB, // Displacement of bodyB
    const size_t& edge_id,                      // In bodyB
    double& toi);

} // namespace ccd
