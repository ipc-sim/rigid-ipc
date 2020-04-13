#pragma once
#include "normal.hpp"

namespace ccd {
namespace geometry {

    template <typename T>
    inline Eigen::Vector2<T> segment_normal(
        const Eigen::Vector2<T>& segment_start,
        const Eigen::Vector2<T>& segment_end,
        bool normalized)
    {
        Eigen::Vector2<T> segment_vec = segment_end - segment_start;
        Eigen::Vector2<T> normal(-segment_vec.y(), segment_vec.x());
        return normalized ? normal.normalized() : normal;
    }

    template <typename T>
    inline Eigen::Vector3<T> triangle_normal(
        const Eigen::Vector3<T>& face_vertex0,
        const Eigen::Vector3<T>& face_vertex1,
        const Eigen::Vector3<T>& face_vertex2,
        bool normalized)
    {
        Eigen::Vector3<T> normal =
            (face_vertex1 - face_vertex0).cross(face_vertex2 - face_vertex0);
        if (normalized) {
            // .normalized() has an inequality check that breaks with Intervals
            normal = normal / normal.norm();
        }
        return normal;
    }

} // namespace geometry
} // namespace ccd
