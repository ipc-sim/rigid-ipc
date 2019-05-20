/** Data structures for impacts between different geometry. */

#ifndef IMPACT_H
#define IMPACT_H

#include <Eigen/Core>
#include <vector>

namespace ccd {

/// Data structure representing an impact of an edge and vertex.
struct EdgeVertexImpact {
    double time;      ///< Time of impact
    int edge_index;   ///< Impacted edge
    double alpha;     ///< Parameter along the edge where the impact occured
    int vertex_index; ///< Impacting vertex

    bool operator==(const EdgeVertexImpact& other) const;
};

/// A vector of edge-vertex impact pointers.
typedef std::vector<EdgeVertexImpact> EdgeVertexImpacts;

/// Data structure representing an impact of an edge and edge.
struct EdgeEdgeImpact {
    /// Time of impact
    double time;
    /// Impacted edge
    int impacted_edge_index;
    /// Parameter along the impacted edge where the impact occured
    double impacted_alpha;
    /// Impacting edge
    int impacting_edge_index;
    /// Parameter along the impacting edge where the impact occured.
    double impacting_alpha;
};

/// A vector of edge-edge impact pointers.
typedef std::vector<EdgeEdgeImpact> EdgeEdgeImpacts;

/**
 * Compare two impacts to determine if impact0 comes before impact1.
 *
 * @param impact0 Impact to check if is came first.
 * @param impact1 Impact to check if is came second.
 * @return A boolean for if impact0.time <= impact1.time.
 */
template <typename Impact>
bool compare_impacts_by_time(Impact impact0, Impact impact1);

/**
 * Convert all edge-vertex impacts to correspoding edge-edge impacts. There may
 * be multiple edge-edge impacts per edge-vertex impact depending on the
 * connectivity.
 *
 * @param edges The matrix of edges where each row is two indices for the
 *     endpoints in the vertices matrix (not a prameter).
 * @param ev_impacts Vector of edge-vertex impacts to convert to edge-edge
 *     impacts.
 * @param ee_impact Vector for where to store the converted edge-edge impacts.
 *     The vector will be cleared.
 * @return Stores the converted edge-edge impacts in ee_impacts.
 */
void convert_edge_vertex_to_edge_edge_impacts(const Eigen::MatrixX2i& edges,
    const EdgeVertexImpacts& ev_impacts, EdgeEdgeImpacts& ee_impacts);

/**
 * Convert all edge-edge impacts to correspoding edge-vertex impacts. There may
 * be multiple edge-vertex impacts per edge-edge impact depending on the
 * connectivity.
 *
 * @param edges A matrix where rows are edges and columns are vertex indices.
 * @param ee_impact Edge-edge impacts to convert to an edge-vertex impacts.
 * @param ev_impact Vector for where to store the converted edge-vertex
 *     impacts. The vector will be cleared.
 * @return Stores the converted edge-vertex impacts in ev_impacts.
 */
void convert_edge_edge_to_edge_vertex_impacts(const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpacts& ee_impacts, EdgeVertexImpacts& ev_impacts);

} // namespace ccd

#endif
