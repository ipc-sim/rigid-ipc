#include "sinc.hpp"

namespace ccd {

double sinc(const double& x)
{
    // We use these epsilons because for example 1 + x^2 = 1 for x < eps_sqrt.
    constexpr static const double eps =
        2.220446049250313080847263336181640625e-16;
    constexpr static const double eps_sqrt = 1.490116119384765625e-8;
    constexpr static const double eps_sqrt_sqrt = 1.220703125e-4;

    if (abs(x) > eps_sqrt_sqrt) {
        return sin(x) / x;
    }

    // approximation by taylor series in x at 0 up to order 1
    double result = 1;

    if (abs(x) >= eps) {
        double x2 = x * x;

        // approximation by taylor series in x at 0 up to order 3
        result -= x2 / 6.0;

        if (abs(x) >= eps_sqrt) {
            // approximation by taylor series in x at 0 up to order 5
            result += (x2 * x2) / 120.0;
        }
    }

    return result;
}

Interval sinc(const Interval& x)
{
    constexpr static const double eps = 1.220703125e-4; // sqrt(sqrt(ϵ))

    // If the domain does not include the challenging interval, compute the
    // value directly.
    if (!overlap(x, Interval(-eps, eps))) {
        return sin(x) / x;
    }

    // Split x in to x_{<ϵ}, x_ϵ, x_{>ϵ}, some might be empty.
    Interval x_lt_eps = Interval(x.lower(), -eps);
    Interval x_gt_eps = Interval(eps, x.upper());

    // Compute x_ϵ and y = sinc(x_{<ϵ}) ∪ sinc(x_{>ϵ})
    Interval y, x_eps;
    if (empty(x_lt_eps) && empty(x_lt_eps)) {
        x_eps = x;
        y = Interval::empty();
    } else if (empty(x_lt_eps)) {
        x_eps = Interval(x.lower(), eps);
        y = sin(x_gt_eps) / x_gt_eps;
    } else if (empty(x_gt_eps)) {
        x_eps = Interval(-eps, x.upper());
        y = sin(x_lt_eps) / x_lt_eps;
    } else {
        // The max range for x_eps is 1
        return hull(
            hull(sin(x_lt_eps) / x_lt_eps, sin(x_gt_eps) / x_gt_eps),
            Interval(1, 1));
    }

    // approximation by taylor series in x at 0 up to order 5
    Interval x_eps2 = x_eps * x_eps;
    return hull(y, x_eps2 * (x_eps2 / 120.0 - 1.0 / 6.0) + 1.0);
}

inline double dsinc_over_x(double x)
{
    static const double eps = 1e-4;

    double x2 = x * x;
    if (abs(x) > eps) {
        return (x * cos(x) - sin(x)) / (x2 * x);
    }

    // approximation by taylor series in x at 0 up to order 5
    return x2 * (-x2 / 840.0 + 1.0 / 30.0) - 1.0 / 3.0;
}

// Compute ∇sinc(||x||) for x ∈ Rⁿ
Eigen::VectorX3d sinc_normx_grad(const Eigen::VectorX3d& x)
{
    return dsinc_over_x(x.norm()) * x;
}

inline double ddsinc_over_x2_minus_dsinc_over_x3(double x)
{
    static const double eps = 0.1;

    double x2 = x * x;
    double x4 = x2 * x2;
    if (abs(x) > eps) {
        return ((3 - x2) * sin(x) + 3 * x * cos(x)) / (x4 * x);
    }

    // approximation by taylor series in x at 0 up to order 5
    return x4 / 7560.0 - x2 / 210.0 + 1.0 / 15.0;
}

// Compute ∇²sinc(||x||) for x ∈ Rⁿ
Eigen::MatrixXX3d sinc_normx_hess(const Eigen::VectorX3d& x)
{
    double normx = x.norm();
    return ddsinc_over_x2_minus_dsinc_over_x3(normx) * x * x.transpose()
        + dsinc_over_x(normx) * Eigen::MatrixXd::Identity(x.size(), x.size());
}

} // namespace ccd
