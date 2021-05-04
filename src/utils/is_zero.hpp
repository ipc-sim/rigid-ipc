#pragma once

#include <interval/interval.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

template <typename T> inline bool is_zero(T x) { return x == 0.0; }

template <> inline bool is_zero(Interval x)
{
    return boost::numeric::zero_in(x);
}

template <int dim, int max_dim>
inline bool is_zero(Vector<Interval, dim, max_dim> X)
{
    // Check that all components contain zero
    // WARNING: This is conservative, but no exact
    for (int i = 0; i < X.size(); i++) {
        if (!is_zero(X(i))) {
            return false;
        }
    }
    return true;
}

} // namespace ipc::rigid
