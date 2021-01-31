#pragma once

#include <vector>

#include <Eigen/Core>

namespace ipc::rigid {

/// @brief Data structure representing an impact of an edge and vertex.
struct EdgeVertexImpact {
    double time;     ///< @brief Time of impact
    long edge_index; ///< @brief Impacted edge
    double alpha; ///< @brief Parameter along the edge where the impact occured
    long vertex_index; ///< @brief Impacting vertex

    EdgeVertexImpact(
        double time, long edge_index, double alpha, long vertex_index);

    bool operator==(const EdgeVertexImpact& other) const;
};

/// @brief Data structure representing an impact of an edge and edge.
struct EdgeEdgeImpact {
    double time;              ///< @brief Time of impact
    long impacted_edge_index; ///< @brief Impacted edge
    double impacted_alpha; ///< @brief Parameter along the impacted edge where
                           ///< the impact occured
    long impacting_edge_index; ///< @brief Impacting edge
    double impacting_alpha; ///< @brief Parameter along the impacting edge where
                            ///< the impact occured.

    EdgeEdgeImpact(
        double time,
        long impacted_edge_index,
        double impacted_alpha,
        long impacting_edge_index,
        double impacting_alpha);

    size_t impacting_node() const { return (impacting_alpha < 0.5 ? 0 : 1); }

    bool operator==(const EdgeEdgeImpact& other) const;
};

/// @brief Data structure representing an impact of an vertex and face.
struct FaceVertexImpact {
    double time;     ///< @brief Time of impact
    long face_index; ///< @brief Impacted face
    double u;        ///< @brief First barycentric coordinate of the face where
                     ///< the impact occured
    double v;        ///< @brief Second barycentric coordinate of the face where
                     ///< the impact occured
    long vertex_index; ///< @brief Impacting vertex

    FaceVertexImpact(
        double time, long face_index, double u, double v, long vertex_index);

    bool operator==(const FaceVertexImpact& other) const;
};

struct Impacts {
    std::vector<EdgeVertexImpact> ev_impacts;
    std::vector<EdgeEdgeImpact> ee_impacts;
    std::vector<FaceVertexImpact> fv_impacts;

    size_t size() const
    {
        return ev_impacts.size() + ee_impacts.size() + fv_impacts.size();
    }

    void clear()
    {
        ev_impacts.clear();
        ee_impacts.clear();
        fv_impacts.clear();
    }
};

/**
 * @brief Compare two impacts to determine if impact0 comes before impact1.
 *
 * @param impact0 Impact to check if is came first.
 * @param impact1 Impact to check if is came second.
 * @return A boolean for if impact0.time <= impact1.time.
 */
template <typename Impact>
bool compare_impacts_by_time(Impact impact1, Impact impact2)
{
    return impact1.time < impact2.time;
}

/**
 * @brief Convert all edge-vertex impacts to correspoding edge-edge impacts.
 *
 * There may be multiple edge-edge impacts per edge-vertex impact depending on
 * the connectivity.
 *
 * @param edges The matrix of edges where each row is two indices for the
 *     endpoints in the vertices matrix (not a prameter).
 * @param ev_impacts Vector of edge-vertex impacts to convert to edge-edge
 *     impacts.
 * @param ee_impact Vector for where to store the converted edge-edge impacts.
 *     The vector will be cleared.
 * @return Stores the converted edge-edge impacts in ee_impacts.
 */
void convert_edge_vertex_to_edge_edge_impacts(
    const Eigen::MatrixX2i& edges,
    const std::vector<EdgeVertexImpact>& ev_impacts,
    std::vector<EdgeEdgeImpact>& ee_impacts);

/**
 * @brief Convert all edge-edge impacts to correspoding edge-vertex impacts.
 *
 * There may be multiple edge-vertex impacts per edge-edge impact depending on
 * the connectivity.
 *
 * @param edges A matrix where rows are edges and columns are vertex indices.
 * @param ee_impact Edge-edge impacts to convert to an edge-vertex impacts.
 * @param ev_impact Vector for where to store the converted edge-vertex
 *     impacts. The vector will be cleared.
 * @return Stores the converted edge-vertex impacts in ev_impacts.
 */
void convert_edge_edge_to_edge_vertex_impacts(
    const Eigen::MatrixX2i& edges,
    const std::vector<EdgeEdgeImpact>& ee_impacts,
    std::vector<EdgeVertexImpact>& ev_impacts);

} // namespace ipc::rigid
