#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc::rigid {

Eigen::SparseMatrix<double> multiply(
    const Eigen::VectorXd& x,
    const std::vector<Eigen::SparseMatrix<double>>& y);

} // namespace ipc::rigid
