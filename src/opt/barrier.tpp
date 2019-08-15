#pragma once
#include "barrier.hpp"

namespace ccd {

namespace opt {

    // Function that grows to infinity as x approaches 0 from the right.
    template <typename T> T spline_barrier(T x, double s)
    {
        if (x <= 0)
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
