#pragma once

#include <Eigen/Core>

namespace ccd {

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_edge_closest_point(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v1);

template <typename T>
T edge_edge_closest_point(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v1);

template <typename T>
T point_triangle_closest_point(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v2);

} // namespace ccd

#include "distance.tpp"
