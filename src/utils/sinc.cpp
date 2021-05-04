#include "sinc.hpp"

namespace ipc::rigid {

// We use these bounds because for example 1 + x^2 = 1 for x < sqrt(ϵ).
static const double taylor_0_bound = std::numeric_limits<double>::epsilon();
static const double taylor_2_bound = sqrt(taylor_0_bound);
static const double taylor_n_bound = sqrt(taylor_2_bound);

double sinc(const double& x)
{
    if (abs(x) >= taylor_n_bound) {
        return sin(x) / x;
    }

    // approximation by taylor series in x at 0 up to order 1
    double result = 1;

    if (abs(x) >= taylor_0_bound) {
        double x2 = x * x;

        // approximation by taylor series in x at 0 up to order 3
        result -= x2 / 6.0;

        if (abs(x) >= taylor_2_bound) {
            // approximation by taylor series in x at 0 up to order 5
            result += (x2 * x2) / 120.0;
        }
    }

    return result;
}

// WARNING: Assumes x is a single value and uses interval arithmetic to account
// for rounding.
Interval _sinc_interval_taylor(double x_double)
{
    Interval x(x_double);

    if (abs(x.lower()) >= taylor_n_bound) {
        return sin(x) / x;
    }
    // approximation by taylor series in x at 0 up to order 5
    // 1 - x^2 / 6 + x^4 / 120 = 1 + x^2 / 6 * (x^2 / 20 - 1)
    Interval x2 = x * x;
    return 1.0 + x2 / 6.0 * (x2 / 20.0 - 1.0);
}

Interval sinc(const Interval& x)
{
    // Define two regions and use even symmetry of sinc.
    // A bound on sinc where it is monotonic ([0, ~4.4934])
    static const double monotonic_bound = 4.4934094579;

    Interval y = Interval::empty(), x_pos = x;
    if (x.lower() < 0) {
        if (x.upper() <= 0) {
            return sinc(-x); // sinc is an even function
        }
        // Split
        y = sinc(Interval(0, -x.lower()));
        x_pos = Interval(0, x.upper());
    }

    // Split the domain into two interval:
    // 1. x ∩ [0, monotonic_bound]
    // 2. x ∩ [monotonic_bound, ∞)

    // Case 1 (Monotonic):
    // WARNING: The following does not account for rounding properly
    Interval x_gt_monotonic = x_pos;
    if (x_pos.lower() <= monotonic_bound) {
        Interval x_monotonic = x_pos;
        if (x_monotonic.upper() > monotonic_bound) {
            x_monotonic = Interval(x_monotonic.lower(), monotonic_bound);
            x_gt_monotonic = Interval(monotonic_bound, x_pos.upper());
        } else {
            x_gt_monotonic = Interval::empty();
        }
        // sinc is monotonically decreasing, so flip a and b.
        // TODO: Set rounding modes here to avoid round off error.
        y = hull(
            y,
            Interval(
                // https://www.wolframalpha.com/input/?i=min+sin%28x%29%2Fx+between+4+and+5
                std::max(
                    _sinc_interval_taylor(x_monotonic.upper()).lower(),
                    -0.271724), // <-- A conservative lower bound for sinc(x)
                std::min(
                    _sinc_interval_taylor(x_monotonic.lower()).upper(), 1.0)));
    }

    // Case 2 (Not necessarily monotonic):
    if (!empty(x_gt_monotonic)) {
        // x_gt_monotonic is larger than one, so the division should be well
        // behaved.
        y = hull(y, sin(x_gt_monotonic) / x_gt_monotonic);
    }

    return y;
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
VectorMax3d sinc_normx_grad(const VectorMax3d& x)
{
    return dsinc_over_x(x.norm()) * x;
}

inline double ddsinc_over_x2_minus_dsinc_over_x3(double x)
{
    static const double eps = 0.1;

    double x2 = x * x;
    double x4 = x2 * x2;
    if (abs(x) > eps) {
        return ((3 - x2) * sin(x) - 3 * x * cos(x)) / (x4 * x);
    }

    // approximation by taylor series in x at 0 up to order 5
    return x4 / 7560.0 - x2 / 210.0 + 1.0 / 15.0;
}

// Compute ∇²sinc(||x||) for x ∈ Rⁿ
MatrixMax3d sinc_normx_hess(const VectorMax3d& x)
{
    double normx = x.norm();
    return ddsinc_over_x2_minus_dsinc_over_x3(normx) * x * x.transpose()
        + dsinc_over_x(normx) * MatrixMax3d::Identity(x.size(), x.size());
}

} // namespace ipc::rigid
