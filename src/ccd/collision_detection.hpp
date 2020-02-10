/**
 * Detect collisions between different geometry.
 * Includes continous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#pragma once

#include <Eigen/Core>
#include <array>

#include <ccd/hash_grid.hpp>
#include <ccd/impact.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {

/// @brief Possible methods for detecting all edge vertex collisions.
enum DetectionMethod {
    BRUTE_FORCE, ///< @brief Use brute-force to detect all collisions.
    HASH_GRID ///< @brief Use a spatial data structure to detect all collisions.
};

namespace CollisionType {
    static const int EDGE_VERTEX = 1;
    static const int EDGE_EDGE = 2;
    static const int FACE_VERTEX = 4;
} // namespace CollisionType

static const char* DetectionMethodNames[] = { "DetectionMethod::HASH_GRID",
                                              "HASH_GRID" };

/**
 * @brief Find all collisions in one time step.
 *
 * @param[in] vertices        The vertices of the bodies.
 * @param[in] displacements   The displacements of the vertices in one
 *                            time-step. There must be an equal number of
 *                            vertices and displacements. The trajectories are
 *                            linear over one time-step and the velocity is
 *                            constant.
 * @param[in] edges           The edges of the bodies defined as pairs of
 *                            indices into the rows of the vertices matrix. Each
 *                            row is an edge.
 * @param[in] faces           The faces of the bodies defined as triplets of
 *                            indices into the rows of the vertices matrix. Each
 *                            row is a face.
 * @param[in] group_ids       If two vertices share a group they are not
 *                            considered possible collisions.
 * @param[in] collision_types Flags for which type of collisions to detect.
 * @param[out] ev_impacts     Store the detected edge-vertex impacts here.
 * @param[out] ee_impacts     Store the detected edge-edge impacts here.
 * @param[out] fv_impacts     Store the detected face-vertex impacts here.
 * @param[in] method          Which method should be used to detect the
 *                            collisions.
 */
void detect_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& faces,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID);

/// @brief Backward compatibility with old definition.
void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID);

/// @brief Backward compatibility with old definition.
void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method = DetectionMethod::HASH_GRID);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_collision_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    const double inflation_radius = 0.0);

/// @brief Backward compatibility with old definition.
void detect_edge_vertex_collision_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    DetectionMethod method = DetectionMethod::HASH_GRID,
    const double inflation_radius = 0.0);

/**
 * @brief Use a brute force method to create a set of all candidate collisions.
 *
 * @param[in] vertices       The vertices of the bodies.
 * @param[in] edges          The edges of the bodies defined as pairs of
 *                           indices into the rows of the vertices matrix. Each
 *                           row is an edge.
 * @param[in] group_ids      If two vertices share a group they are not
 *                           considered possible collisions.
 * @param[out] ev_candidates Vector of candidates to build.
 * @param[out] ee_candidates Vector of candidates to build.
 * @param[out] fv_candidates Vector of candidates to build.
 */
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates);

/**
 * @brief Use a hash grid method to create a set of all candidate collisions.
 *
 * @param[in] vertices       The vertices of the bodies.
 * @param[in] displacements  The displacements of the vertices in one time-step.
 *                           There must be an equal number of vertices and
 *                           displacements. The trajectories are linear over one
 *                           time-step and the velocity is constant.
 * @param[in] edges          The edges of the bodies defined as pairs of indices
 *                           into the rows of the vertices matrix. Each row is
 *                           an edge.
 * @param[in] group_ids      If two vertices share a group they are not
 *                           considered possible collisions.
 * @param[out] ev_candidates Vector of candidates to build.
 */
void detect_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const EdgeVertexCandidates& ev_candidates,
    const EdgeEdgeCandidates& ee_candidates,
    const FaceVertexCandidates& fv_candidates,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts);

/**
 * @brief Determine if a single edge-vertext pair intersects.
 *
 * If the edge and vertex are impacting, store a new EdgeVertexImpact in
 * ev_impacts.
 *
 * @param[in] Vi           The position, in 2D, of the first endpoint of the
 *                         edge.
 * @param[in] Vj           The position, in 2D, of the second endpoint of the
 *                         edge.
 * @param[in] Vk           The position, in 2D, of the vertex.
 * @param[in] Ui           The displacement of \f$V_i(t)\f$ over the time-step.
 * @param[in] Uj           The displacement of \f$V_j(t)\f$ over the time-step.
 * @param[in] Uk           The displacement of \f$V_k(t)\f$ over the time-step.
 * @param[in] edge_id      Index of the edge.
 * @param[in] vertex_id    Index of the vertex.
 * @param[out] ev_impacts  List of impacts on to which new impacts are pushed.
 */
void detect_edge_vertex_collisions_narrow_phase(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const EdgeVertexCandidate& ev_candidate,
    EdgeVertexImpacts& ev_impacts);

void detect_edge_edge_collisions_narrow_phase(
    const Eigen::VectorXd& Vi,
    const Eigen::VectorXd& Vj,
    const Eigen::VectorXd& Vk,
    const Eigen::VectorXd& Vl,
    const Eigen::VectorXd& Ui,
    const Eigen::VectorXd& Uj,
    const Eigen::VectorXd& Uk,
    const Eigen::VectorXd& Ul,
    const EdgeEdgeCandidate& ee_candidate,
    EdgeEdgeImpacts& ee_impacts);

void detect_face_vertex_collisions_narrow_phase(
    const Eigen::VectorXd& Vi,
    const Eigen::VectorXd& Vj,
    const Eigen::VectorXd& Vk,
    const Eigen::VectorXd& Vl,
    const Eigen::VectorXd& Ui,
    const Eigen::VectorXd& Uj,
    const Eigen::VectorXd& Uk,
    const Eigen::VectorXd& Ul,
    const FaceVertexCandidate& fv_candidate,
    FaceVertexImpacts& fv_impacts);

} // namespace ccd
