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
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& edge_id,                 // In bodyB
    double& toi);

/// Find time-of-impact between two rigid bodies
bool compute_edge_edge_time_of_impact(
    const physics::RigidBody& bodyA,
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& edgeA_id,                // In bodyA
    const physics::RigidBody& bodyB,
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& edgeB_id,                // In bodyB
    double& toi);

/// Find time-of-impact between two rigid bodies
bool compute_face_vertex_time_of_impact(
    const physics::RigidBody& bodyA,
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& face_id,                 // In bodyB
    double& toi);

} // namespace ccd
