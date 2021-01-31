#pragma once

#include <Eigen/Core>

namespace ipc::rigid {

template <typename T>
inline void flatten(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& x)
{
    x.resize(x.size(), 1);
}

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> flat(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& x)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vx = x;
    vx.resize(vx.size(), 1);
    return vx;
}

template <typename T>
inline void unflatten(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& x, const int cols)
{
    assert(x.size() % cols == 0);
    x.resize(x.size() / cols, cols);
}
} // namespace ipc::rigid
