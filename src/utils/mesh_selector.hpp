#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <vector>

namespace ipc::rigid {

class MeshSelector {
public:
    MeshSelector() {}
    MeshSelector(
        size_t num_vertices,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F);

    size_t vertex_to_edge(size_t vi) const { return m_vertex_to_edge[vi]; }
    size_t vertex_to_face(size_t vi) const { return m_vertex_to_face[vi]; }

    size_t face_to_edge(size_t fi, size_t fi_ei) const
    {
        return m_faces_to_edges(fi, fi_ei);
    }

    const Eigen::MatrixXi& face_to_edges() const { return m_faces_to_edges; }

    size_t edge_to_face(size_t ei) const { return m_edge_to_face[ei]; }

    size_t codim_vertices_to_vertices(size_t vi) const
    {
        return m_codim_vertices_to_vertices[vi];
    }

    size_t codim_edges_to_edges(size_t ei) const
    {
        return m_codim_edges_to_edges[ei];
    }

    const std::vector<size_t>& codim_edges_to_edges() const
    {
        return m_codim_edges_to_edges;
    }

    size_t num_codim_vertices() const
    {
        return m_codim_vertices_to_vertices.size();
    }

    size_t num_codim_edges() const { return m_codim_edges_to_edges.size(); }

protected:
    /// Map from vertices to the minimum index edge containing the vertex.
    void init_vertex_to_edge(size_t num_vertices, const Eigen::MatrixXi& E);

    /// Map from vertices to the minimum index face containing the vertex.
    void init_vertex_to_face(size_t num_vertices, const Eigen::MatrixXi& F);

    /// Each row is a face with the three indices of the edges.
    void init_face_to_edge(
        const Eigen::MatrixXi& F,
        const Eigen::SparseMatrix<size_t>& vertices_to_edge);

    void init_edge_to_face(const size_t num_edges);

    void init_codim_vertices_to_vertices(
        size_t num_vertices, const Eigen::MatrixXi& F);

    void init_codim_edges_to_edges(size_t num_edges);

    /// A mapping from a vertex id to the smallest adjacent face id.
    std::vector<size_t> m_vertex_to_edge;
    /// A mapping from a vertex id to the smallest adjacent face id.
    std::vector<size_t> m_vertex_to_face;
    /// A mapping from a edge id to the smallest adjacent face id.
    std::vector<size_t> m_edge_to_face;
    /// A mapping from a face to the indices of the face's edges.
    Eigen::MatrixXi m_faces_to_edges;
    /// A map from codimensional vertices ids to full vertices ids.
    std::vector<size_t> m_codim_vertices_to_vertices;
    /// A map from codimensional edges ids to full edges ids.
    std::vector<size_t> m_codim_edges_to_edges;
};

} // namespace ipc::rigid
