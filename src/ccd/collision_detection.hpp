/**
 * Detect collisions between different geometry.
 * Includes continous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#pragma once

#include <array>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <ccd/hash_grid.hpp>
#include <ccd/impact.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {

/// @brief Possible methods for detecting all edge vertex collisions.
enum DetectionMethod {
    BRUTE_FORCE, ///< @brief Use brute-force to detect all collisions.
    HASH_GRID ///< @brief Use a spatial data structure to detect all collisions.
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    DetectionMethod,
    { { HASH_GRID, "hash_grid" }, { BRUTE_FORCE, "brute_force" } });

namespace CollisionType {
    static const int EDGE_VERTEX = 1;
    static const int EDGE_EDGE = 2;
    static const int FACE_VERTEX = 4;
} // namespace CollisionType

/**
 * @brief Find all collisions in one time step.
 *
 * @param[in] vertices_t0     The vertices of the bodies at the start of the
 *                            time-step.
 * @param[in] vertices_t1     The vertices of the bodies at the end of the
 *                            time-step.
 * @param[in] edges           The edges of the bodies defined as pairs of
 *                            indices into the rows of the vertices matrix. Each
 *                            row is an edge.
 * @param[in] faces           The faces of the bodies defined as triplets of
 *                            indices into the rows of the vertices matrix. Each
 *                            row is a face.
 * @param[in] group_ids       If two vertices share a group they are not
 *                            considered possible collisions.
 * @param[in] collision_types Flags for which type of collisions to detect.
 * @param[out] impacts        Store the detected impacts here.
 * @param[in] method          Which method should be used to detect the
 *                            collisions.
 */
void detect_collisions(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID);

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
    Candidates& candidates,
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

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Candidates& candidates,
    ConcurrentImpacts& impacts);

/// @brief Determine if an edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const Eigen::Vector2d& edge_vertex0_t0,
    const Eigen::Vector2d& edge_vertex1_t0,
    const Eigen::Vector2d& vertex_t0,
    const Eigen::Vector2d& edge_vertex0_t1,
    const Eigen::Vector2d& edge_vertex1_t1,
    const Eigen::Vector2d& vertex_t1,
    double& toi,
    double& alpha);

bool detect_edge_edge_collisions_narrow_phase(
    const Eigen::Vector3d& edge0_vertex0_t0,
    const Eigen::Vector3d& edge0_vertex1_t0,
    const Eigen::Vector3d& edge1_vertex0_t0,
    const Eigen::Vector3d& edge1_vertex1_t0,
    const Eigen::Vector3d& edge0_vertex0_t1,
    const Eigen::Vector3d& edge0_vertex1_t1,
    const Eigen::Vector3d& edge1_vertex0_t1,
    const Eigen::Vector3d& edge1_vertex1_t1,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha);

bool detect_face_vertex_collisions_narrow_phase(
    const Eigen::Vector3d& face_vertex0_t0,
    const Eigen::Vector3d& face_vertex1_t0,
    const Eigen::Vector3d& face_vertex2_t0,
    const Eigen::Vector3d& vertex_t0,
    const Eigen::Vector3d& face_vertex0_t1,
    const Eigen::Vector3d& face_vertex1_t1,
    const Eigen::Vector3d& face_vertex2_t1,
    const Eigen::Vector3d& vertex_t1,
    double& toi,
    double& u,
    double& v);

} // namespace ccd
