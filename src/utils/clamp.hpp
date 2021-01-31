#pragma once

#include <interval/interval.hpp>

namespace ipc::rigid {

template <typename T> inline T clamp_to_01(const T& x)
{
    return x > T(1.0) ? T(1.0) : (x < T(0.0) ? T(0.0) : x);
}

template <> inline Interval clamp_to_01(const Interval& x)
{
    return boost::numeric::intersect(x, Interval(0, 1));
}

} // namespace ipc::rigid
