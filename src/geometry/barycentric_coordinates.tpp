// Copied this function from libigl and templated it.
#pragma once
#include "barycentric_coordinates.hpp"

#include <logger.hpp>

namespace ccd {
namespace geometry {

    // Computes the barycentric coordinates of a point in a triangle.
    //
    // Computes the barycentric coordinates of any point in the plane
    // containing the triangle.
    template <typename T>
    inline void barycentric_coordinates(
        const Eigen::Vector3<T>& point,
        const Eigen::Vector3<T>& triangle_vertex0,
        const Eigen::Vector3<T>& triangle_vertex1,
        const Eigen::Vector3<T>& triangle_vertex2,
        T& u,
        T& v,
        T& w)
    {
        const Eigen::Array<T, 3, 1> v0 =
            triangle_vertex1.array() - triangle_vertex0.array();
        const Eigen::Array<T, 3, 1> v1 =
            triangle_vertex2.array() - triangle_vertex0.array();
        const Eigen::Array<T, 3, 1> v2 =
            point.array() - triangle_vertex0.array();

        T d00 = (v0 * v0).sum();
        T d01 = (v0 * v1).sum();
        T d11 = (v1 * v1).sum();
        T d20 = (v2 * v0).sum();
        T d21 = (v2 * v1).sum();

        T denom = d00 * d11 - d01 * d01;
        if (denom == 0) {
            spdlog::error("Unable to compute barycentric coordinates!");
            assert(denom != 0);
        }

        v = (d11 * d20 - d01 * d21) / denom;
        w = (d00 * d21 - d01 * d20) / denom;
        u = 1.0 - (v + w);
    }

} // namespace geometry
} // namespace ccd
