#pragma once

#include <Eigen/Core>

namespace ccd {

template <typename T>
inline void flatten(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& x)
{
    x.resize(x.size(), 1);
}

template <typename T>
inline void unflatten(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& x, const int cols)
{
    assert(x.size() % cols == 0);
    x.resize(x.size() / cols, cols);
}
} // namespace ccd
