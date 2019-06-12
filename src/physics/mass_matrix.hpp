#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace physics {

    void mass_matrix(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges, Eigen::SparseMatrix<double>& mass_matrix);

    void mass_vector(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges, Eigen::VectorXd& mass_vector);
}
}
