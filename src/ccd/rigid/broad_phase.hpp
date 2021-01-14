#pragma once

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <interval/interval.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Discrete Collision Detection
// NOTE: Yes, I know this is inside the CCD directory.
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Continous Collision Detection
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/// @brief Use a hash grid method to create a set of all candidate collisions.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

template <typename T>
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<T>& poses,
    const int collision_types,
    const size_t bodyA_id,
    const size_t bodyB_id,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

typedef tbb::enumerable_thread_specific<ipc::Candidates>
    ThreadSpecificCandidates;

void merge_local_candidates(
    const ThreadSpecificCandidates& storages, ipc::Candidates& candidates);

inline ipc::AABB vertex_aabb(
    const Eigen::MatrixXd& V, const size_t vi, double inflation_radius = 0)
{
    return ipc::AABB(
        V.row(vi).array() - inflation_radius,
        V.row(vi).array() + inflation_radius);
}

inline ipc::AABB vertex_aabb(
    const Eigen::MatrixXI& V, const size_t vi, double inflation_radius = 0)
{
    const auto& v = V.row(vi);
    Eigen::VectorX3d min(v.size());
    Eigen::VectorX3d max(v.size());
    for (int i = 0; i < v.size(); i++) {
        min(i) = (v(i) - inflation_radius).lower();
        max(i) = (v(i) + inflation_radius).upper();
    }
    return ipc::AABB(min, max);
}

template <typename T>
inline ipc::AABB edge_aabb(
    const Eigen::MatrixX<T>& V,
    const Eigen::Vector2i& E,
    double inflation_radius = 0)
{
    return ipc::AABB(
        vertex_aabb(V, E(0), inflation_radius),
        vertex_aabb(V, E(1), inflation_radius));
}

template <typename T>
inline ipc::AABB face_aabb(
    const Eigen::MatrixX<T>& V,
    const Eigen::Vector3i& F,
    double inflation_radius = 0)
{
    return ipc::AABB(
        vertex_aabb(V, F(0), inflation_radius),
        ipc::AABB(
            vertex_aabb(V, F(1), inflation_radius),
            vertex_aabb(V, F(2), inflation_radius)));
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Intersection Detection
///////////////////////////////////////////////////////////////////////////////

/// Use a BVH to create a set of all candidate intersections.
void detect_intersection_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    std::vector<ipc::EdgeFaceCandidate>& ef_candidates);

} // namespace ccd
