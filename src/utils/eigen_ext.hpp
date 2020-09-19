#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ipc/utils/eigen_ext.hpp>

namespace Eigen {

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x);

template <typename T> inline Matrix2<T> Hat(T x);
template <typename T> inline Matrix3<T> Hat(Vector3<T> x);
template <typename T> inline MatrixXX3<T> Hat(VectorX3<T> x);

} // namespace Eigen

#include "eigen_ext.tpp"
