// Copied this function from libigl and templated it.
#pragma once

#include <Eigen/Core>

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace geometry {

    /// Computes the barycentric coordinates of a point in a triangle.
    ///
    /// Computes the barycentric coordinates of any point in the plane
    /// containing the triangle.
    template <typename T>
    inline void barycentric_coordinates(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2,
        T& u,
        T& v,
        T& w);

} // namespace geometry
} // namespace ccd

#include "barycentric_coordinates.tpp"
