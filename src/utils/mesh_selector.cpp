#include "mesh_selector.hpp"

#include <Eigen/SparseCore>
#include <igl/Timer.h>

namespace ipc::rigid {

Eigen::SparseMatrix<size_t>
map_vertex_pairs_to_edge(size_t num_vertices, const Eigen::MatrixXi& E)
{
    // TODO: this could be faster using a non std unordered map
    std::vector<Eigen::Triplet<size_t>> vertices_to_edge_triplets;
    vertices_to_edge_triplets.reserve(E.rows());

    for (size_t i = 0; i < E.rows(); i++) {
        vertices_to_edge_triplets.emplace_back(
            std::min(E(i, 0), E(i, 1)), std::max(E(i, 0), E(i, 1)), i);
    }

    Eigen::SparseMatrix<size_t> vertices_to_edge(num_vertices, num_vertices);
    vertices_to_edge.setFromTriplets(
        vertices_to_edge_triplets.begin(), vertices_to_edge_triplets.end());
    vertices_to_edge.makeCompressed();

    return vertices_to_edge;
}

MeshSelector::MeshSelector(
    size_t num_vertices, const Eigen::MatrixXi& E, const Eigen::MatrixXi& F)
{
    init_vertex_to_edge(num_vertices, E);
    init_vertex_to_face(num_vertices, F);

    if (F.rows()) {
        auto vertices_to_edge = map_vertex_pairs_to_edge(num_vertices, E);
        init_face_to_edge(F, vertices_to_edge);
        init_edge_to_face(E.rows());
    }

    init_codim_vertices_to_vertices(num_vertices, E);
    init_codim_edges_to_edges(E.rows());
}

// Map from vertices to the minimum index edge containing the vertex.
void MeshSelector::init_vertex_to_edge(
    size_t num_vertices, const Eigen::MatrixXi& E)
{
    m_vertex_to_edge.resize(num_vertices, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < E.rows(); i++) {
        for (size_t j = 0; j < E.cols(); j++) {
            m_vertex_to_edge[E(i, j)] = std::min(m_vertex_to_edge[E(i, j)], i);
        }
    }
}

// Map from vertices to the minimum index face containing the vertex.
void MeshSelector::init_vertex_to_face(
    size_t num_vertices, const Eigen::MatrixXi& F)
{
    m_vertex_to_face.resize(num_vertices, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < F.rows(); i++) {
        for (size_t j = 0; j < F.cols(); j++) {
            m_vertex_to_face[F(i, j)] = std::min(m_vertex_to_face[F(i, j)], i);
        }
    }
}

// Each row is a face with the three indices of the edges.
void MeshSelector::init_face_to_edge(
    const Eigen::MatrixXi& F,
    const Eigen::SparseMatrix<size_t>& vertices_to_edge)
{
    m_faces_to_edges.resize(F.rows(), 3);
    for (size_t i = 0; i < F.rows(); i++) {
        for (size_t j = 0; j < F.cols(); j++) {
            m_faces_to_edges(i, j) = vertices_to_edge.coeff(
                std::min(F(i, j), F(i, (j + 1) % 3)),
                std::max(F(i, j), F(i, (j + 1) % 3)));
        }
    }
}

void MeshSelector::init_edge_to_face(const size_t num_edges)
{
    assert(m_faces_to_edges.size());

    if (!num_edges) {
        return;
    }

    m_edge_to_face.resize(num_edges, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < m_faces_to_edges.rows(); i++) {
        for (size_t j = 0; j < m_faces_to_edges.cols(); j++) {
            m_edge_to_face[m_faces_to_edges(i, j)] =
                std::min(m_edge_to_face[m_faces_to_edges(i, j)], i);
        }
    }
}

void MeshSelector::init_codim_vertices_to_vertices(
    size_t num_vertices, const Eigen::MatrixXi& E)
{
    std::vector<bool> is_vertex_codim(num_vertices, true);
    for (size_t i = 0; i < E.rows(); i++) {
        for (size_t j = 0; j < E.cols(); j++) {
            is_vertex_codim[E(i, j)] = false;
        }
    }
    size_t num_codim_vertices =
        std::count(is_vertex_codim.begin(), is_vertex_codim.end(), true);
    m_codim_vertices_to_vertices.reserve(num_codim_vertices);
    for (size_t i = 0; i < is_vertex_codim.size(); i++) {
        if (is_vertex_codim[i]) {
            m_codim_vertices_to_vertices.push_back(i);
        }
    }
}

void MeshSelector::init_codim_edges_to_edges(size_t num_edges)
{
    std::vector<bool> is_edge_codim(num_edges, true);
    for (size_t i = 0; i < m_faces_to_edges.rows(); i++) {
        for (size_t j = 0; j < m_faces_to_edges.cols(); j++) {
            is_edge_codim[m_faces_to_edges(i, j)] = false;
        }
    }
    size_t num_codim_edges =
        std::count(is_edge_codim.begin(), is_edge_codim.end(), true);
    m_codim_edges_to_edges.reserve(num_codim_edges);
    for (size_t i = 0; i < is_edge_codim.size(); i++) {
        if (is_edge_codim[i]) {
            m_codim_edges_to_edges.push_back(i);
        }
    }
}

} // namespace ipc::rigid
