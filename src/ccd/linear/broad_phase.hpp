/**
 * Detect collisions between different geometry.
 * Includes continuous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/broad_phase/hash_grid.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates_linear(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/**
 * @brief Use a brute force method to create a set of all candidate collisions.
 *
 * @param[in] vertices     The vertices of the bodies.
 * @param[in] edges        The edges of the bodies defined as pairs of indices
 *                         into the rows of the vertices matrix. Each row is an
 *                         edge.
 * @param[in] group_ids    If two vertices share a group they are not considered
 *                         possible collisions.
 * @param[out] candidates  Vectors of candidates to build.
 */
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates);

/**
 * @brief Use a hash grid method to create a set of all candidate collisions.
 *
 * @param[in] vertices_t0     The vertices of the bodies at the start of the
 *                            time-step.
 * @param[in] vertices_t1     The vertices of the bodies at the end of the
 *                            time-step.
 * @param[in] edges          The edges of the bodies defined as pairs of indices
 *                           into the rows of the vertices matrix. Each row is
 *                           an edge.
 * @param[in] group_ids      If two vertices share a group they are not
 *                           considered possible collisions.
 * @param[out] candidates    Vector of candidates to build.
 */
void detect_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

/// @brief Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_linear_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

inline AABB vertex_aabb(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const size_t vi,
    double inflation_radius)
{
    return AABB(
        V0.row(vi).cwiseMin(V1.row(vi)).array() - inflation_radius,
        V0.row(vi).cwiseMax(V1.row(vi)).array() + inflation_radius);
}

inline AABB edge_aabb(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::Vector2i& E,
    double inflation_radius)
{
    return AABB(
        vertex_aabb(V0, V1, E(0), inflation_radius),
        vertex_aabb(V0, V1, E(1), inflation_radius));
}

inline AABB face_aabb(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::Vector3i& F,
    double inflation_radius)
{
    return AABB(
        vertex_aabb(V0, V1, F(0), inflation_radius),
        AABB(
            vertex_aabb(V0, V1, F(1), inflation_radius),
            vertex_aabb(V0, V1, F(2), inflation_radius)));
}

void detect_collision_candidates_linear_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    const size_t bodyA_id,
    const size_t bodyB_id,
    Candidates& candidates,
    const double inflation_radius);

} // namespace ipc::rigid
