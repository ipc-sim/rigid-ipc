#pragma once
#include "barrier.hpp"
#include <autodiff/autodiff.h>
#include <utils/not_implemented_error.hpp>
namespace ccd {

namespace opt {

    template <typename T> T ipc_barrier(T d, double dhat)
    {
        if (d <= T(0))
            return T(std::numeric_limits<double>::infinity());
        if (d >= dhat)
            return T(0);
        // b(d) = -(d-d̂)²ln(d / d̂)
        T dhat_T = T(dhat);
        return -(d - dhat_T) * (d - dhat_T) * log(d / dhat_T);
    }

    template <typename T> T poly_log_barrier(T x, double s)
    {
        if (x <= T(0))
            return T(std::numeric_limits<double>::infinity());
        if (x >= s)
            return T(0);

        // y:= x/eps;  -log(y)*(2*y^3-3*y^2+1)
        T y = x / s;
        T g = -log(y) * ((2 * y - 3) * y * y + 1);
        return g;
    }

    template <typename T> T spline_barrier(T x, double s)
    {
        if (x <= T(0))
            return T(std::numeric_limits<double>::infinity());
        if (x >= s)
            return T(0);
        T x_s = x / s;
        // g(x) = (x / s)^3 - 3 * (x / s)^2 + 3 * (x / s)
        T g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        return 1 / g - 1;
    }

} // namespace opt
} // namespace ccd
