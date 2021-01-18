#include "rigid_trajectory_aabb.hpp"

namespace ccd {

typedef physics::Pose<Interval> PoseI;

Eigen::VectorX3I vertex_trajectory_aabb(
    const physics::RigidBody& body,
    const PoseI& pose_t0, // Pose of body at t=0
    const PoseI& pose_t1, // Pose of body at t=1
    size_t vertex_id,     // In body
    const Interval& t)
{
    // Compute the pose at time t
    PoseI pose = PoseI::interpolate(pose_t0, pose_t1, t);
    // Get the world vertex of the edges at time t
    return body.world_vertex(pose, vertex_id);
}

Eigen::VectorX3I edge_trajectory_aabb(
    const physics::RigidBody& body,
    const PoseI& pose_t0, // Pose of body at t=0
    const PoseI& pose_t1, // Pose of body at t=1
    size_t edge_id,       // In body
    const Interval& t,
    const Interval& alpha)
{
    // Compute the pose at time t
    PoseI pose = PoseI::interpolate(pose_t0, pose_t1, t);
    // Get the world vertex of the edges at time t
    Eigen::VectorX3I e0 = body.world_vertex(pose, body.edges(edge_id, 0));
    Eigen::VectorX3I e1 = body.world_vertex(pose, body.edges(edge_id, 1));
    return (e1 - e0) * alpha + e0;
}

Eigen::VectorX3I face_trajectory_aabb(
    const physics::RigidBody& body,
    const PoseI& pose_t0, // Pose of body at t=0
    const PoseI& pose_t1, // Pose of body at t=1
    size_t face_id,       // In body
    const Interval& t,
    const Interval& u,
    const Interval& v)
{
    // Compute the pose at time t
    PoseI pose = PoseI::interpolate(pose_t0, pose_t1, t);
    // Get the world vertex of the edges at time t
    Eigen::VectorX3I f0 = body.world_vertex(pose, body.faces(face_id, 0));
    Eigen::VectorX3I f1 = body.world_vertex(pose, body.faces(face_id, 1));
    Eigen::VectorX3I f2 = body.world_vertex(pose, body.faces(face_id, 2));
    return (f1 - f0) * u + (f2 - f0) * v + f0;
}

Eigen::VectorX3I edge_vertex_aabb(
    const physics::RigidBody& bodyA, // Body of the vertex
    const PoseI& poseA_t0,           // Pose of bodyA at t=0
    const PoseI& poseA_t1,           // Pose of bodyA at t=1
    size_t vertex_id,                // In bodyA
    const physics::RigidBody& bodyB, // Body of the edge
    const PoseI& poseB_t0,           // Pose of bodyB at t=0
    const PoseI& poseB_t1,           // Pose of bodyB at t=1
    size_t edge_id,                  // In bodyB
    const Interval& t,
    const Interval& alpha)
{
    return vertex_trajectory_aabb(bodyA, poseA_t0, poseA_t1, vertex_id, t)
        - edge_trajectory_aabb(bodyB, poseB_t0, poseB_t1, edge_id, t, alpha);
}

Eigen::VectorX3I edge_edge_aabb(
    const physics::RigidBody& bodyA, // Body of the first edge
    const PoseI& poseA_t0,           // Pose of bodyA at t=0
    const PoseI& poseA_t1,           // Pose of bodyA at t=1
    size_t edgeA_id,                 // In bodyA
    const physics::RigidBody& bodyB, // Body of the second edge
    const PoseI& poseB_t0,           // Pose of bodyB at t=0
    const PoseI& poseB_t1,           // Pose of bodyB at t=1
    size_t edgeB_id,                 // In bodyB
    const Interval& t,
    const Interval& alpha,
    const Interval& beta)
{
    return edge_trajectory_aabb(bodyA, poseA_t0, poseA_t1, edgeA_id, t, alpha)
        - edge_trajectory_aabb(bodyB, poseB_t0, poseB_t1, edgeB_id, t, beta);
}

Eigen::VectorX3I face_vertex_aabb(
    const physics::RigidBody& bodyA,         // Body of the vertex
    const physics::Pose<Interval>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<Interval>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,                        // In bodyA
    const physics::RigidBody& bodyB,         // Body of the triangle
    const physics::Pose<Interval>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<Interval>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id,                          // In bodyB
    const Interval& t,
    const Interval& u,
    const Interval& v)
{
    return vertex_trajectory_aabb(bodyA, poseA_t0, poseA_t1, vertex_id, t)
        - face_trajectory_aabb(bodyB, poseB_t0, poseB_t1, face_id, t, u, v);
}

} // namespace ccd
