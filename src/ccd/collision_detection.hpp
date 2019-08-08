/**
 * Detect collisions between different geometry.
 * Includes continous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

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

static const char* DetectionMethodNames[] = { "BRUTE_FORCE", "HASH_GRID" };

/**
 * @brief Find all edge-vertex collisions in one time step.
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
 * @param[out] ev_impacts    Reference for the vector for where to store the
 *                           detected impacts.
 * @param[in] method         Which method should be used to detect the
 *                           collisions.
 * @param[in] reset_impacts  Reset the vector of impacts.
 */
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method = BRUTE_FORCE,
    bool reset_impacts = true);

/// @brief Backward compatibility with old definition.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method = BRUTE_FORCE,
    bool reset_impacts = true);

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

/// @brief Use broad-phase method to create a set of candidate collisions.
void detect_edge_vertex_collision_candidates(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    DetectionMethod method,
    const double inflation_radius = 0.0);

/**
 * @brief Use a brute force method to create a set of all candidate edge-vertex
 * collisions.
 *
 * @param[in] vertices       The vertices of the bodies.
 * @param[in] displacements  The displacements of the vertices in one
 *                           time-step. There must be an equal number of
 *                           vertices and displacements. The trajectories are
 *                           linear over one time-step and the velocity is
 *                           constant.
 * @param[in] edges          The edges of the bodies defined as pairs of
 *                           indices into the rows of the vertices matrix. Each
 *                           row is an edge.
 * @param[in] group_ids      If two vertices share a group they are not
 *                           considered possible collisions.
 * @param[out] ev_candidates Vector of candidates to build.
 */
void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    const double inflation_radius = 0.0);

/**
 * @brief Use a hash grid method to create a set of all candidate edge-vertex
 * collisions.
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
void detect_edge_vertex_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    const double inflation_radius = 0.0);

///////////////////////////////////////////////////////////////////////////////
// Narrow-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_edge_vertex_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const EdgeVertexCandidates& ev_candidates,
    EdgeVertexImpacts& ev_impacts);

void detect_edge_vertex_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const EdgeVertexCandidates& ev_candidates,
    const Eigen::MatrixXb& skip_pair,
    EdgeVertexImpacts& ev_impacts);

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
void detect_edge_vertex_collisions_narrow_phase(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const EdgeVertexCandidate& ev_candidate,
    EdgeVertexImpacts& ev_impacts);

/**
 * @brief Compute the time of impact of a point and edge moving in 2D.
 *
 * @param[in] Vi      The position, in 2D, of the first endpoint of the edge.
 * @param[in] Vj      The position, in 2D, of the second endpoint of the edge.
 * @param[in] Vk      The position, in 2D, of the vertex.
 * @param[in] Ui      The displacement of \f$V_i(t)\f$ over the time-step.
 * @param[in] Uj      The displacement of \f$V_j(t)\f$ over the time-step.
 * @param[in] Uk      The displacement of \f$V_k(t)\f$ over the time-step.
 * @param[out] toi    The time of impact, \f$\tau_I \in [0, 1]\f$, where the
 *                    start of the time-step is at \f$t=0\f$ and the end of the
 *                    time-step is at time \f$t=1\f$. The computed output value
 *                    is stored in this variable.
 * @param[out] alpha  The spatial position of the impact along the edge,
 *                    \f$\alpha \in [0, 1]\f$, where the start of the edge is at
 *                    \f$\alpha=0\f$ and the end of the edge is at
 *                    \f$\alpha=1\f$. The computed output value is stored in
 *                    this variable.
 *
 * @return  A boolean for if a collision occured in time 0 to 1. If there is a
 *          collision the computed time of impact is stored in toi and the
 *          impact points spatial parameter along the edge is stored in alpha.
 */
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    double& toi,
    double& alpha,
    const double tolerance = 1e-8);

/**
 * @brief Convert a temporal parameterization to a spatial parameterization.
 *
 * Move the vertices to the time of impact and then computing the ratio
 * \f$\frac{V_k - V_i}{V_j - V_i}\f$. This function is best defined for a time
 * of impact where vertex0 lies on the line spaned by the edge vertices.
 *
 * @param[in] Vi      The first end-point of the edge.
 * @param[in] Vj      The second end-point of the edge.
 * @param[in] Vk      The vertex to find the spatial parameterization along the
 *                    edge at time t.
 * @param[in] Ui      The displacement of \f$V_i(t)\f$ over the time-step.
 * @param[in] Uj      The displacement of \f$V_j(t)\f$ over the time-step.
 * @param[in] Uk      The displacement of \f$V_k(t)\f$ over the time-step.
 * @param[in] t       The time of impact to move the vertices.
 *
 * @param[out] alpha  The spatial parameterization, \f$s\f$, such that
 *                    \f$V_k(t) = (V_j(t) - V_i(t))*s + V_i(t)\f$. The computed
 *                     output value is stored in this variable.
 *
 * @return A boolean for if their is a valid spatial parameterization. Stores
 *         the valid spatial parameter in alpha.
 */
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const double t,
    double& alpha,
    const double tolerance = 1e-8);

} // namespace ccd
#endif
