#pragma once

#include <interval/interval.hpp>
#include <physics/rigid_body.hpp>

namespace ipc::rigid {

VectorMax3I vertex_trajectory_aabb(
    const RigidBody& body,
    const Pose<Interval>& pose_t0, // Pose of body at t=0
    const Pose<Interval>& pose_t1, // Pose of body at t=1
    size_t vertex_id,              // In body
    const Interval& t = Interval(0, 1));

VectorMax3I edge_trajectory_aabb(
    const RigidBody& body,
    const Pose<Interval>& pose_t0, // Pose of body at t=0
    const Pose<Interval>& pose_t1, // Pose of body at t=1
    size_t edge_id,                // In body
    const Interval& t = Interval(0, 1),
    const Interval& alpha = Interval(0, 1));

VectorMax3I face_trajectory_aabb(
    const RigidBody& body,
    const Pose<Interval>& pose_t0, // Pose of body at t=0
    const Pose<Interval>& pose_t1, // Pose of body at t=1
    size_t face_id,                // In body
    const Interval& t = Interval(0, 1),
    const Interval& u = Interval(0, 1),
    const Interval& v = Interval(0, 1));

VectorMax3I edge_vertex_aabb(
    const RigidBody& bodyA,         // Body of the vertex
    const Pose<Interval>& poseA_t0, // Pose of bodyA at t=0
    const Pose<Interval>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,               // In bodyA
    const RigidBody& bodyB,         // Body of the edge
    const Pose<Interval>& poseB_t0, // Pose of bodyB at t=0
    const Pose<Interval>& poseB_t1, // Pose of bodyB at t=1
    size_t edge_id,                 // In bodyB
    const Interval& t = Interval(0, 1),
    const Interval& alpha = Interval(0, 1));

VectorMax3I edge_edge_aabb(
    const RigidBody& bodyA,         // Body of the first edge
    const Pose<Interval>& poseA_t0, // Pose of bodyA at t=0
    const Pose<Interval>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,                // In bodyA
    const RigidBody& bodyB,         // Body of the second edge
    const Pose<Interval>& poseB_t0, // Pose of bodyB at t=0
    const Pose<Interval>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id,                // In bodyB
    const Interval& t = Interval(0, 1),
    const Interval& alpha = Interval(0, 1),
    const Interval& beta = Interval(0, 1));

VectorMax3I face_vertex_aabb(
    const RigidBody& bodyA,         // Body of the vertex
    const Pose<Interval>& poseA_t0, // Pose of bodyA at t=0
    const Pose<Interval>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,               // In bodyA
    const RigidBody& bodyB,         // Body of the triangle
    const Pose<Interval>& poseB_t0, // Pose of bodyB at t=0
    const Pose<Interval>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id,                 // In bodyB
    const Interval& t = Interval(0, 1),
    const Interval& u = Interval(0, 1),
    const Interval& v = Interval(0, 1));

} // namespace ipc::rigid
