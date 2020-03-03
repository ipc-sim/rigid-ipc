#pragma once

#include <Eigen/Core>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace geometry {

    template <typename T>
    inline Eigen::Vector2<T> segment_normal(
        const Eigen::Vector2<T>& segment_start,
        const Eigen::Vector2<T>& segment_end,
        bool normalized = true);

    template <typename T>
    inline Eigen::Vector3<T> triangle_normal(
        const Eigen::Vector3<T>& face_vertex0,
        const Eigen::Vector3<T>& face_vertex1,
        const Eigen::Vector3<T>& face_vertex2,
        bool normalized = true);

} // namespace geometry
} // namespace ccd

#include "normal.tpp"
