// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequlity constraints on a function.
#include "barrier.hpp"

#include <barrier/barrier_chorner.hpp>
#include <logger.hpp>

namespace ccd {

namespace opt {

    // template spetialization
    template <> double spline_barrier<double>(double x, double s)
    {
        return barrier_horner_compensated(x, s);
    }

    double barrier_gradient(double x, double s, BarrierType barrier_type)
    {
        switch (barrier_type) {
        case BarrierType::IPC:
            return ipc_barrier_gradient(x, s);
        case BarrierType::POLY_LOG:
            return poly_log_barrier_gradient(x, s);
        case BarrierType::SPLINE:
            return spline_barrier_gradient(x, s);
        }
    }

    double ipc_barrier_gradient(double d, double dhat)
    {
        if (d <= 0.0 || d >= dhat) {
            return 0.0;
        }
        // b(d) = -(d - d̂)²ln(d / d̂)
        // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
        //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
        return (d - dhat) * (-2 * log(d / dhat) - (d - dhat) / d);
    }

    double poly_log_barrier_gradient(double x, double s)
    {
        // (6x/s² - 6x²/s³)log(x/s) + (-1 + 3x²/s² - 2x³/s³)/x
        if (x <= 0)
            return 0;

        if (x >= s)
            return 0.0;

        double x2 = x * x, s2 = s * s;
        double x3 = x2 * x, s3 = s2 * s;

        return (6 * x / s2 - 6 * x2 / s3) * log(x / s)
            + (-1 + 3 * x2 / s2 - 2 * x3 / s3) / x;
    }

    // Derivative of the spline_barrier function with respect to x.
    double spline_barrier_gradient(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        double x_s = x / s;
        // g(x) = (x / s)³ - 3 * (x / s)² + 3 * (x / s)
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        // g'(x) = 3 * x² / s³ - 6 * x / s² + 3 / s
        //       = (3 * (x / s)² - 6 * (x / s) + 3 ) / s
        double dg = (3 + x_s * (-6 + 3 * x_s)) / s; // Horner's method
        // ϕ'(x) = -g⁻² * g'
        return -1 / (g * g) * dg;
    }

    // Second derivative of the spline_barrier function with respect to x.
    double spline_barrier_hessian(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        double x_s = x / s;
        // g(x) = (x / s)³ - 3 * (x / s)² + 3 * (x / s)
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        // g'(x) = 3 * x² / s³ - 6 * x / s² + 3 / s
        //       = (3 * (x / s)² - 6 * (x / s) + 3 ) / s
        double dg = (3 + x_s * (-6 + 3 * x_s)) / s; // Horner's method
        // g''(x) = 6 * x / s³ - 6 / s² = (6 * (x / s) - 6) / s²
        double ddg = (6 * x_s - 6) / (s * s);
        // ϕ''(x) = 2g^-3 * (g')² + -g⁻² * g'' = (2 / g * (g')² - g'') / g²
        return (2 / g * dg * dg - ddg) / (g * g);
    }

} // namespace opt
} // namespace ccd
