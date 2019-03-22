/**
 * Detect collisions between different geometry.
 * Includes continous collision detection to compute the time of impact.
 * Supported geometry: point vs edge
 */

#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include <Eigen/Core>
#include <array>

#include <ccd/impact.hpp>

namespace ccd {

/// Possible methods for detecting all edge vertex collisions.
enum DetectionMethod {
    BRUTE_FORCE, ///< Use brute-force to detect all collisions.
    HASH_MAP     ///< Use a spatial data structure to detect all collisions.
};

static const char* DetectionMethodNames[] = { "BRUTE_FORCE", "HASH_MAP" };

/**
 * Convert a temporal parameterization to a spatial parameterization.
 *
 * Move the vertices to the time of impact and then computing the ratio
 * \f$\frac{p_0 - p_1}{p_2 - p_1}\f$. This function is best defined for a time
 * of impact where vertex0 lies on the line spaned by the edge vertices.
 *
 * @param[in] vertex0 The vertex to find the spatial parameterization along the
 *                    edge at time t.
 * @param[in] displacment0 The displacement of vertex0 over the time-step (i.e.
 *                         vertex0' = vertex0 + t * displacement0).
 * @param[in] edge_vertex1 The first end-point of the edge.
 * @param[in] edge_displacement1 The displacement of edge_vertex1 over the
 *                               time-step.
 * @param[in] edge_vertex2 The second end-point of the edge.
 * @param[in] edge_displacement2 The displacement of edge_vertex2 over the
 *                               time-step.
 * @param[in] t The time of impact to move the vertices.
 * @param[out] alpha The spatial parameterization, \f$s\f$, such that
 *                   \f$p_0(t) = (p_2(t) - p_1(t))*s + p_1(t)\f$. The computed
 *                   output value is stored in this variable.
 * @return A boolean for if their is a valid spatial parameterization. Stores
 *         the valid spatial parameter in alpha.
 */
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& vertex0,
    const Eigen::Vector2d& displacement0, const Eigen::Vector2d& edge_vertex1,
    const Eigen::Vector2d& edge_displacement1,
    const Eigen::Vector2d& edge_vertex2,
    const Eigen::Vector2d& edge_displacement2, const double t, double& alpha);

/**
 * Compute the time of impact of a point and edge moving in 2D.
 *
 * @param[in] vertex0 The position, in 2D, of the vertex.
 * @param[in] displacment0 The displacement, in 2D, of the vertex.
 * @param[in] edge_vertex1 The position, in 2D, of the first endpoint of the
 * edge.
 * @param[in] edge_displacement1 The displacement, in 2D, of the first endpoint
 * of the edge.
 * @param[in] edge_vertex1 The position, in 2D, of the second endpoint of the
 * edge.
 * @param[in] edge_displacement1 The displacement, in 2D, of the second endpoint
 * of the edge.
 * @param[out] toi The time of impact, \f$\tau_I \in [0, 1]\f$, where the start
 * of the time-step is at \f$t=0\f$ and the end of the time-step is at time
 * \f$t=1\f$. The computed output value is stored in this variable.
 * @param[out] alpha The spatial position of the impact along the edge,
 * \f$\alpha \in [0, 1]\f$, where the start of the edge is at \f$\alpha=0\f$ and
 * the end of the edge is at time \f$\alpha=1\f$. The computed output value is
 * stored in this variable.
 * @return A boolean for if a collision occured in time 0 to 1. If there is a
 * collision the computed time of impact is stored in toi and the impact points
 * spatial parameter along the edge is stored in alpha.
 */
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& vertex0,
    const Eigen::Vector2d& displacement0, const Eigen::Vector2d& edge_vertex1,
    const Eigen::Vector2d& edge_displacement1,
    const Eigen::Vector2d& edge_vertex2,
    const Eigen::Vector2d& edge_displacement2, double& toi, double& alpha);

/**
 * Find all edge-vertex collisions in one time step.
 *
 * @param[in] vertices The vertices of the bodies.
 * @param[in] displacements The displacements of the vertices in one time-step.
 * There must be an equal number of vertices and displacments. The trajectories
 * are linear over one time-step and the velocity is constant.
 * @param[in] edges The edges of the bodies defined as pairs of indices into the
 * rows of the vertices matrix. Each row is an edge.
 * @param[out] ev_impacts Reference for the vector for where to store the
 * detected impacts.
 * @param[in] method Which method should be used to detect the collisions.
 * @return All impacts as impact structures containing the impacting edge and
 * vertex index and the time of impact are stored in ev_impacts.
 */
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts, DetectionMethod method = BRUTE_FORCE);

/**
 * Find all edge-vertex collisions in one time step using brute-force
 * comparisons of all edges and all vertices.
 *
 * @param[in] vertices The vertices of the bodies.
 * @param[in] displacements The displacements of the vertices in one time-step.
 * There must be an equal number of vertices and displacments. The trajectories
 * are linear over one time-step and the velocity is constant.
 * @param[in] edges The edges of the bodies defined as pairs of indices into the
 * rows of the vertices matrix. Each row is an edge.
 * @param[out] ev_impacts Reference for the vector for where to store the
 * detected impacts.
 * @return All impacts as impact structures containing the impacting edge and
 * vertex index and the time of impact are stored in ev_impacts.
 */
void detect_edge_vertex_collisions_brute_force(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts);

/**
 * Find all edge-vertex collisions in one time step using spatial-hashing to
 * only compare points and edge in the same cells.
 *
 * TODO: Remove the [[noreturn]]  after implementing this function
 * The [[noreturn]] is used to supress warnings.
 *
 * @param[in] vertices The vertices of the bodies.
 * @param[in] displacements The displacements of the vertices in one time-step.
 * There must be an equal number of vertices and displacments. The trajectories
 * are linear over one time-step and the velocity is constant.
 * @param[in] edges The edges of the bodies defined as pairs of indices into the
 * rows of the vertices matrix. Each row is an edge.
 * @param[out] ev_impacts Reference for the vector for where to store the
 * detected impacts.
 * @return All impacts as impact structures containing the impacting edge and
 * vertex index and the time of impact are stored in ev_impacts.
 */
[[noreturn]] void detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges, EdgeVertexImpacts& ev_impacts);

} // namespace ccd
#endif
