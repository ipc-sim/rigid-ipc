#include "finitediff.hpp"

#include <array>
#include <vector>

/**
  Based on the functions on https://github.com/PatWie/CppNumericalSolvers
  and rewritten to use Eigen
  **/

namespace ccd {
bool compare_gradient(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{

    assert(x.rows() == y.rows());

    for (long d = 0; d < x.rows(); ++d) {
        double scale = std::max(std::max(fabs(x[d]), fabs(y[d])), double(1.0));
        if (fabs(x[d] - y[d]) > 1e-6 * scale)
            return false;
    }
    return true;
}

void finite_gradient(
    const Eigen::VectorXd& x, std::function<double(const Eigen::VectorXd&)> value,
    Eigen::VectorXd& grad, size_t accuracy)
{
    // accuracy can be 0, 1, 2, 3
    assert(accuracy <= 3);

    const double eps = 2.2204e-6;
    // clang-format off
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    // clang-format on
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
}
