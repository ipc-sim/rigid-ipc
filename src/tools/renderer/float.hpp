#pragma once

#include <Eigen/Core>

namespace swr {

using Float = double;

typedef Eigen::Matrix<Float, 3, 3> Matrix3F;
typedef Eigen::Matrix<Float, 4, 4> Matrix4F;
typedef Eigen::Matrix<Float, 2, 1> Vector2F;
typedef Eigen::Matrix<Float, 3, 1> Vector3F;
typedef Eigen::Matrix<Float, 4, 1> Vector4F;

} // namespace swr
