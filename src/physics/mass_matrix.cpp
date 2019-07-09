#include "mass_matrix.hpp"

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace physics {

    // Construct the sparse mass matrix for the given mesh (V, E).
    void mass_matrix(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        Eigen::SparseMatrix<double>& mass_matrix)
    {
        Eigen::VectorXd vertex_masses;
        mass_vector(vertices, edges, vertex_masses);
        mass_matrix = Eigen::SparseDiagonal<double>(vertex_masses);
    }

    // Construct the diagonal of the mass matrix for the given mesh (V, E).
    void mass_vector(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        Eigen::VectorXd& vertex_masses)
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

    double total_mass(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
    {
        double mass = 0.0;
        for (long i = 0; i < edges.rows(); i++) {
            double edge_length
                = (vertices.row(edges(i, 1)) - vertices.row(edges(i, 0)))
                      .norm();
            mass += edge_length;
        }
        return mass;
    }
} // namespace physics
} // namespace ccd
