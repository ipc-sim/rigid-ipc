#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace tensor {

    Eigen::SparseMatrix<double> multiply(const Eigen::VectorXd& x,
        const std::vector<Eigen::SparseMatrix<double>>& y);
} // namespace tensor
} // namespace ccd
