#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace physics {

    /// @brief Construct the sparse mass matrix for the given mesh (V, E).
    void mass_matrix(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        Eigen::SparseMatrix<double>& mass_matrix);

    /// @brief Construct the diagonal of the mass matrix for the given mesh
    /// (V, E).
    void mass_vector(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        Eigen::VectorXd& mass_vector);

    /// @brief Construct the total mass for the given mesh
    double total_mass(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges);
} // namespace physics
} // namespace ccd
