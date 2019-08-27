/**
  Based on the functions on https://github.com/PatWie/CppNumericalSolvers
  and rewritten to use Eigen
  **/
#ifndef FINITEDIFF_H
#define FINITEDIFF_H

#include <Eigen/Core>

namespace ccd {

/**
Enumeration of available orders of accuracy for finite differences. The
corresponding integer values are used internally and should be ignored.
*/
enum AccuracyOrder {
    SECOND = 0, /*!< Second order accuracy. */
    FOURTH = 1, /*!< Fourth order accuracy. */
    SIXTH = 2,  /*!< Sixth order accuracy. */
    EIGHTH = 3  /*!< Eighth order accuracy. */
};

/**
Compare if two gradients are close enough.

@param x The first gradient to compare.
@param y The second gradient to compare against.
@return A boolean for if x and y are close to the same value.
*/
bool compare_gradient(const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const double test_eps = 1e-4,
    const std::string& msg = "compare_gradient ");
bool compare_jacobian(const Eigen::MatrixXd& x,
    const Eigen::MatrixXd& y,
    const double test_eps = 1e-4,
    const std::string& msg = "compare_jacobian ");
/**
Compute the gradient of a function at a point using finite differences.

@param x The points in \f$\mathbb{R}^{n \times m}\f$ at which to compute the
    gradient.
@param value The function that goes from \f$\mathbb{R}^{n \times m} \rightarrow
    \mathbb{R}\f$
@param grad Where to store the gradient of the function.
@param accuracy The accuracy of the finite differences (Options: 0, 1, 2, 3)
@return The computed gradient is stored in the grad parameter.
*/
void finite_gradient(const Eigen::VectorXd& x,
    std::function<double(const Eigen::VectorXd&)> value,
    Eigen::VectorXd& grad,
    AccuracyOrder accuracy = SECOND,
    const double eps = 1.0e-8);

void finite_jacobian(const Eigen::VectorXd& x,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> value,
    Eigen::MatrixXd& jac,
    const AccuracyOrder accuracy = SECOND,
    const double eps = 1.0e-8);
} // namespace ccd
#endif
