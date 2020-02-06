#pragma once

#include <Eigen/Core>

namespace ccd {

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_edge_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v1);

template <typename T>
T edge_edge_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v1);

template <typename T>
T point_triangle_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v2);

//-----------------------------------------------------------------------------
// Signed Distances
//-----------------------------------------------------------------------------

template <typename T>
T point_edge_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e_v1);

template <typename T>
T edge_edge_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e0_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& e1_v1);

template <typename T>
T point_triangle_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v0,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v1,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& t_v2,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& v);

template <typename T>
T point_plane_signed_distance(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,  ///< Vertex
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& p0, ///< Point on the plane.
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& n   ///< Normal of the plane.
);

} // namespace ccd

#include "distance.tpp"
