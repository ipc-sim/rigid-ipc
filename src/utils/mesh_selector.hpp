#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

namespace ccd {

class MeshSelector {
public:
    MeshSelector() {}
    MeshSelector(
        size_t num_vertices,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F);

    size_t vertex_to_face(size_t vi) const { return m_vertex_to_face[vi]; }
    size_t face_to_edge(size_t fi, size_t fi_ei) const
    {
        return m_face_to_edges(fi, fi_ei);
    }
    const Eigen::MatrixXi& face_to_edges() const { return m_face_to_edges; }
    size_t edge_to_face(size_t ei) const { return m_edge_to_face[ei]; }

protected:
    /// A mapping from a vertex id to the smallest adjacent face id.
    std::vector<size_t> m_vertex_to_edge;
    /// A mapping from a vertex id to the smallest adjacent face id.
    std::vector<size_t> m_vertex_to_face;
    /// A mapping from a edge id to the smallest adjacent face id.
    std::vector<size_t> m_edge_to_face;
    /// A mapping from a face to the indices of the face's edges.
    Eigen::MatrixXi m_face_to_edges;
};

} // namespace ccd
