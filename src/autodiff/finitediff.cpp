#include "finitediff.hpp"

#include <array>
#include <vector>

#include <logger.hpp>

// Based on the functions on https://github.com/PatWie/CppNumericalSolvers
// and rewritten to use Eigen

namespace ccd {

// Compare if two gradients are close enough.
bool compare_gradient(const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const double test_eps,
    const std::string& msg)
{
    assert(x.rows() == y.rows());

    bool same = true;
    for (long d = 0; d < x.rows(); ++d) {
        double scale = std::max(std::max(fabs(x[d]), fabs(y[d])), double(1.0));
        double abs_diff = fabs(x[d] - y[d]);

        if (abs_diff > test_eps * scale) {
            spdlog::warn(
                "{} eps={:.3e} r={} x={:.3e} y={:.3e} |x-y|={:.3e} |x-y|/|x|={:.3e} |x-y|/|y={:3e}",
                msg, test_eps, d, x(d), y(d), abs_diff, abs_diff / fabs(x(d)),
                abs_diff / fabs(y(d)));
            same = false;
        }
    }
    return same;
}

bool compare_jacobian(const Eigen::MatrixXd& x,
    const Eigen::MatrixXd& y,
    const double test_eps,
    const std::string& msg)
{
    assert(x.rows() == y.rows());
    assert(x.cols() == y.cols());

    bool same = true;
    for (long d = 0; d < x.rows(); ++d) {
        for (long c = 0; c < x.cols(); ++c) {
            double scale
                = std::max(std::max(fabs(x(d, c)), fabs(y(d, c))), double(1.0));

            double abs_diff = fabs(x(d, c) - y(d, c));

            if (abs_diff > test_eps * scale) {
                spdlog::warn(
                    "{} eps={:.3e} r={} c={} x={:.3e} y={:.3e} |x-y|={:.3e} |x-y|/|x|={:.3e} |x-y|/|y={:3e}",
                    msg, test_eps, d, c, x(d, c), y(d, c), abs_diff,
                    abs_diff / fabs(x(d, c)), abs_diff / fabs(y(d, c)));
                same = false;
            }
        }
    }
    return same;
}

// Compute the gradient of a function at a point using finite differences.
void finite_gradient(const Eigen::VectorXd& x,
    std::function<double(const Eigen::VectorXd&)> value,
    Eigen::VectorXd& grad,
    AccuracyOrder accuracy,
    const double eps)
{
    // Create an array of the coefficients for finite differences.
    // See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    // clang-format off
    // The external coefficients, c1, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    // The internal coefficients, c2, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    // clang-format on
    // The denominators of the finite difference.
    static const std::array<double, 4> dd = { { 2, 12, 60, 840 } };

    grad.resize(x.rows());

    const size_t innerSteps = 2 * (accuracy + 1);
    const double ddVal = dd[accuracy] * eps;

    Eigen::VectorXd xx = x;
    for (long d = 0; d < x.rows(); d++) {
        grad[d] = 0;
        for (size_t s = 0; s < innerSteps; ++s) {
            double tmp = x[d];
            xx[d] += coeff2[accuracy][s] * eps;
            grad[d] += coeff[accuracy][s] * value(xx);
            xx[d] = tmp;
        }
        grad[d] /= ddVal;
    }
}

void finite_jacobian(const Eigen::VectorXd& x,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> value,
    Eigen::MatrixXd& jac,
    const AccuracyOrder accuracy,
    const double eps)
{
    // Create an array of the coefficients for finite differences.
    // See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    // clang-format off
    // The external coefficients, c1, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    // The internal coefficients, c2, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    // clang-format on
    // The denominators of the finite difference.
    static const std::array<double, 4> dd = { { 2, 12, 60, 840 } };

    jac.resize(value(x).rows(), x.rows());

    const size_t innerSteps = 2 * (accuracy + 1);
    const double ddVal = dd[accuracy] * eps;

    Eigen::VectorXd xx = x;
    for (long d = 0; d < x.rows(); d++) {
        jac.col(d).setZero();
        for (size_t s = 0; s < innerSteps; ++s) {
            double tmp = x[d];
            xx[d] += coeff2[accuracy][s] * eps;
            jac.col(d) += coeff[accuracy][s] * value(xx);
            xx[d] = tmp;
        }
        jac.col(d) /= ddVal;
    }
}
} // namespace ccd
