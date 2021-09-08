#pragma once

#include <ipc/barrier/barrier.hpp>

#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

template <typename T> T barrier(const T& x, double s, BarrierType barrier_type)
{
    switch (barrier_type) {
    case BarrierType::IPC:
        return ipc::barrier(x, s);
    case BarrierType::POLY_LOG:
        return poly_log_barrier(x, s);
    case BarrierType::SPLINE:
        return spline_barrier(x, s);
    default:
        throw NotImplementedError(
            fmt::format("Invalid barrier type: {:d}", int(barrier_type))
                .c_str());
    }
}

template <typename T> T poly_log_barrier(const T& x, double s)
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

template <typename T> T spline_barrier(const T& x, double s)
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

} // namespace ipc::rigid
