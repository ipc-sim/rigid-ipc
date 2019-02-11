/**
  Based on the functions on https://github.com/PatWie/CppNumericalSolvers
  and rewritten to use Eigen
  **/
#ifndef FINITEDIFF_H
#define FINITEDIFF_H

#include <Eigen/Core>

namespace ccd {
bool compare_gradient(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

void finite_gradient(
    const Eigen::VectorXd& x, std::function<double(const Eigen::VectorXd&)> value,
    Eigen::VectorXd& grad, size_t accuracy = 0);
}
#endif
