#pragma once

#include "eigen_ext.hpp"

namespace Eigen {

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x)
{
    return MatrixX<T>(x.asDiagonal()).sparseView();
}

template <typename T> Matrix2<T> Hat(T x)
{
    Matrix2<T> M;
    M.row(0) << T(0), -x;
    M.row(1) << x, T(0);
    return M;
}

template <typename T> Matrix3<T> Hat(Eigen::Vector3<T> x)
{
    Matrix3<T> M;
    M.row(0) << T(0), -x.z(), x.y();
    M.row(1) << x.z(), T(0), -x.x();
    M.row(2) << -x.y(), x.x(), T(0);
    return M;
}

template <typename T> MatrixXX3<T> Hat(Eigen::VectorX3<T> x)
{
    if (x.size() == 1) {
        return Hat(x(0));
    } else {
        return Hat(Eigen::Vector3<T>(x));
    }
}

} // namespace Eigen
