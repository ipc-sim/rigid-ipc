#include "mesh_selector.hpp"

#include <iostream>
#include <unordered_map>

#include <fmt/format.h>

namespace ccd {

MeshSelector::MeshSelector(
    size_t num_vertices, const Eigen::MatrixXi& E, const Eigen::MatrixXi& F)
{
    auto edge_hash = [&num_vertices](const Eigen::Vector2i& e) {
        size_t vi_min = e.minCoeff(), vi_max = e.maxCoeff();
        return vi_min * num_vertices + vi_max;
    };
    auto edges_equal = [](const Eigen::Vector2i& e0,
                          const Eigen::Vector2i& e1) {
        size_t e0_vi_min = e0.minCoeff(), e0_vi_max = e0.maxCoeff();
        size_t e1_vi_min = e1.minCoeff(), e1_vi_max = e1.maxCoeff();
        return e0_vi_min == e1_vi_min && e0_vi_max == e1_vi_max;
    };
    std::unordered_map<
        Eigen::Vector2i, size_t, decltype(edge_hash), decltype(edges_equal)>
        vertices_to_edge(/*bucket_count=*/E.rows(), edge_hash, edges_equal);
    for (size_t i = 0; i < E.rows(); i++) {
        vertices_to_edge.emplace(E.row(i), i);
    }

    // Map from vertices to the minimum index edge containing the vertex.
    m_vertex_to_edge.resize(num_vertices, std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < E.rows(); i++) {
        for (size_t j = 0; j < E.cols(); j++) {
            m_vertex_to_edge[E(i, j)] = std::min(m_vertex_to_edge[E(i, j)], i);
        }
    }

    // Map from vertices to the minimum index face containing the vertex.
    m_vertex_to_face.resize(num_vertices, std::numeric_limits<size_t>::max());
    // Each row is a face with the three indices of the edges.
    m_face_to_edges.resize(F.rows(), 3);
    for (size_t i = 0; i < F.rows(); i++) {
        for (size_t j = 0; j < F.cols(); j++) {
            m_vertex_to_face[F(i, j)] = std::min(m_vertex_to_face[F(i, j)], i);
            m_face_to_edges(i, j) =
                vertices_to_edge[Eigen::Vector2i(F(i, j), F(i, (j + 1) % 3))];
        }
    }

    m_edge_to_face.resize(E.rows(), std::numeric_limits<size_t>::max());
    if (E.rows()) {
        for (size_t i = 0; i < m_face_to_edges.rows(); i++) {
            for (size_t j = 0; j < m_face_to_edges.cols(); j++) {
                m_edge_to_face[m_face_to_edges(i, j)] =
                    std::min(m_edge_to_face[m_face_to_edges(i, j)], i);
            }
        }
    }
}

} // namespace ccd
