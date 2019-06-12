#include "mass_matrix.hpp"

namespace ccd {
namespace physics {

    void mass_matrix(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges, Eigen::SparseMatrix<double>& mass_matrix)
    {
        Eigen::VectorXd vertex_masses;
        mass_vector(vertices, edges, vertex_masses);
        mass_matrix = Eigen::MatrixXd(vertex_masses.asDiagonal()).sparseView();
    }

    void mass_vector(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges, Eigen::VectorXd& vertex_masses)
    {
        vertex_masses = Eigen::VectorXd::Zero(vertices.rows());
        for (long i = 0; i < edges.rows(); i++) {
            double edge_length
                = (vertices.row(edges(i, 1)) - vertices.row(edges(i, 0)))
                      .norm();
            // Add vornoi areas to the vertex weight
            vertex_masses(edges(i, 0)) += edge_length / 2;
            vertex_masses(edges(i, 1)) += edge_length / 2;
        }
    }
} // namespace physics
} // namespace ccd
