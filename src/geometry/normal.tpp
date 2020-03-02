#pragma once
#include "normal.hpp"

namespace ccd {
namespace geometry {

    template <typename T>
    inline Eigen::VectorX3<T> segment_normal(
        const Eigen::VectorX3<T>& segment_start,
        const Eigen::VectorX3<T>& segment_end,
        bool normalized)
    {
        assert(segment_start.size() == 2);
        assert(segment_end.size() == 2);
        Eigen::Vector2<T> segment_vec = segment_end - segment_start;
        Eigen::Vector2<T> normal(-segment_vec.y(), segment_vec.x());
        return normalized ? normal.normalized() : normal;
    }

    template <typename T>
    inline Eigen::VectorX3<T> triangle_normal(
        const Eigen::VectorX3<T>& face_vertex0,
        const Eigen::VectorX3<T>& face_vertex1,
        const Eigen::VectorX3<T>& face_vertex2,
        bool normalized)
    {
        assert(face_vertex0.size() == 3);
        assert(face_vertex1.size() == 3);
        assert(face_vertex2.size() == 3);
        Eigen::Vector3<T> normal =
            ((face_vertex1 - face_vertex0).template head<3>())
                .cross((face_vertex2 - face_vertex0).template head<3>());
        return normalized ? normal.normalized() : normal;
    }

} // namespace geometry
} // namespace ccd
