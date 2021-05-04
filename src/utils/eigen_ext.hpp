#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

typedef Eigen::DiagonalMatrix<double, 3> DiagonalMatrix3d;

template <typename T>
Eigen::SparseMatrix<T> SparseDiagonal(const VectorX<T>& x);

template <typename T> inline Matrix2<T> Hat(T x);
template <typename T> inline Matrix3<T> Hat(Vector3<T> x);
template <typename T> inline MatrixMax3<T> Hat(VectorMax3<T> x);

} // namespace ipc

#include "eigen_ext.tpp"
