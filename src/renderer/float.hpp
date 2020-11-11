#pragma once

#include <Eigen/Core>

namespace swr {

using Float = double;

} // namespace swr

namespace Eigen {
typedef Matrix<swr::Float, 3, 3> Matrix3F;
typedef Matrix<swr::Float, 4, 4> Matrix4F;
typedef Matrix<swr::Float, 2, 1> Vector2F;
typedef Matrix<swr::Float, 3, 1> Vector3F;
typedef Matrix<swr::Float, 4, 1> Vector4F;
} // namespace Eigen
