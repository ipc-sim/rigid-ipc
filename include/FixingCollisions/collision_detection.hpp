/**
    Detection collisions between different geometry.
    Includes continous collision detection to compute the time of impact.
    Supported geometry: point vs edge
*/

#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include <Eigen/Core>
#include <vector>

namespace ccd {

/** Structure representing an impact of an edge and vertex. */
struct Impact {
    /** Impacting vertex. */
    int vertex_index;
    /** Impacted edge. */
    int edge_index;
    /** Time of impact. */
    double time;
};

typedef std::shared_ptr<Impact> ImpactPtr;
typedef std::vector<ImpactPtr> Impacts;
typedef std::shared_ptr<Impacts> ImpactsPtr;

/** Possible methods for detecting all edge vertex collisions. */
enum DetectionMethod { BRUTE_FORCE, HASH_MAP };

/**
    Convert a temporal parameterization to a spatial parameterization.
    This is done by moving the vertices to the time of impact and then
    computing the ratio (p0 - p1) / (p2 - p1). This function is best defined
    for a time of impact where vertex0 lies on the line spaned by the edge
    vertices.

    @param vertex0 The vertex to find the spatial parameterization along the
                   edge at time t.
    @param displacment0 The displacement of vertex0 over the time-step
                        (i.e. vertex0' = vertex0 + t * displacement0).
    @param edge_vertex1 The first end-point of the edge.
    @param edge_displacement1 The displacement of edge_vertex1 over the
                              time-step.
    @param edge_vertex2 The second end-point of the edge.
    @param edge_displacement2 The displacement of edge_vertex2 over the
                              time-step.
    @param t The time of impact to move the vertices.
    @return The spatial parameterization such that
            p0(t) = (p2(t) - p1(t))*s + p1(t).
*/
double temporal_parameterization_to_spatial(const Eigen::VectorXd& vertex0,
    const Eigen::VectorXd& displacement0, const Eigen::VectorXd& edge_vertex1,
    const Eigen::VectorXd& edge_displacement1,
    const Eigen::VectorXd& edge_vertex2,
    const Eigen::VectorXd& edge_displacement2, double t);

/**
    Compute the time of impact of a point and edge moving in 2D.

    @param vertex0 The position, in 2D, of the vertex.
    @param displacment0 The displacement, in 2D, of the vertex.
    @param edge_vertex1 The position, in 2D, of the first endpoint of the
   edge.
    @param edge_displacement1 The displacement, in 2D, of the first endpoint
   of the edge.
    @param edge_vertex1 The position, in 2D, of the second endpoint of the
   edge.
    @param edge_displacement1 The displacement, in 2D, of the second
   endpoint of the edge.
    @return The time of impact, tI in [0, 1] where the start of the time
   step is at time t=0 and the end of the time step is at time t=1.
*/
double compute_edge_vertex_time_of_impact(const Eigen::Vector2d& vertex0,
    const Eigen::Vector2d& displacement0, const Eigen::Vector2d& edge_vertex1,
    const Eigen::Vector2d& edge_displacement1,
    const Eigen::Vector2d& edge_vertex2,
    const Eigen::Vector2d& edge_displacement2);

/**
    Find all edge-vertex collisions in one time step.

    @param vertices The vertices of the bodies.
    @param displacements The displacements of the vertices in one time-step.
                         There must be an equal number of vertices and
                         displacments. The trajectories are linear over one
                         time-step and the velocity is constant.
    @param edges The edges of the bodies defined as pairs of indices into the
                 rows of the vertices matrix. Each row is an edge.
    @param method Which method should be used to detect the collisions.
    @return All impacts as impact structures containing the impacting edge and
            vertex index and the time of impact.
*/
ImpactsPtr detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    DetectionMethod method = BRUTE_FORCE);

/**
    Find all edge-vertex collisions in one time step using brute-force
    comparisons of all edges and all vertices.

    @param vertices The vertices of the bodies.
    @param displacements The displacements of the vertices in one time-step.
                         There must be an equal number of vertices and
                         displacments. The trajectories are linear over one
                         time-step and the velocity is constant.
    @param edges The edges of the bodies defined as pairs of indices into the
                 rows of the vertices matrix. Each row is an edge.
    @return All impacts as impact structures containing the impacting edge and
            vertex index and the time of impact.
*/
ImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges);

/**
    Find all edge-vertex collisions in one time step using spatial-hashing to
    only compare points and edge in the same cells.

    @param vertices The vertices of the bodies.
    @param displacements The displacements of the vertices in one time-step.
                         There must be an equal number of vertices and
                         displacments. The trajectories are linear over one
                         time-step and the velocity is constant.
    @param edges The edges of the bodies defined as pairs of indices into the
                 rows of the vertices matrix. Each row is an edge.
    @return All impacts as impact structures containing the impacting edge and
            vertex index and the time of impact.
*/
ImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges);

} // namespace ccd
#endif
