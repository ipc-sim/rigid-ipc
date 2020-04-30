#pragma once

#include "eigen_ext.hpp"

namespace Eigen {

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x)
{
    return MatrixX<T>(x.asDiagonal()).sparseView();
}

template <typename T> inline Matrix2<T> Hat(T x)
{
    Matrix2<T> M;
    T zero(0);
    M.row(0) << zero, -x;
    M.row(1) << x, zero;
    return M;
}

template <typename T> inline Matrix3<T> Hat(Eigen::Vector3<T> x)
{
    Matrix3<T> M;
    T zero(0);
    M.row(0) << zero, -x.z(), x.y();
    M.row(1) << x.z(), zero, -x.x();
    M.row(2) << -x.y(), x.x(), zero;
    return M;
}

template <typename T> inline MatrixXX3<T> Hat(Eigen::VectorX3<T> x)
{
    if (x.size() == 1) {
        return Hat(x(0));
    } else {
        return Hat(Eigen::Vector3<T>(x));
    }
}

} // namespace Eigen
