/**
 * Detect collisions between different geometry.
 * Includes continuous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>

#include <ccd/ccd.hpp>
#include <ccd/impact.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
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
    ipc::Candidates& candidates);

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
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

} // namespace ccd
