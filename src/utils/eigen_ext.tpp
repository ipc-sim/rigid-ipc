#pragma once

#include "eigen_ext.hpp"

namespace Eigen {

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x)
{
    return MatrixX<T>(x.asDiagonal()).sparseView();
}

template <typename T> MatrixXX3<T> Hat(Eigen::VectorX3<T> x)
{
    assert(x.size() != 2);
    MatrixXX3<T> M;
    if (x.size() == 1) {
        M = Matrix2<T>();
        M.row(0) << 0, -x(0);
        M.row(1) << x(0), 0;
    } else {
        M = Matrix3<T>();
        M.row(0) << 0, -x.z(), x.y();
        M.row(1) << x.z(), 0, -x.x();
        M.row(2) << -x.y(), x.x(), 0;
    }
    return M;
}

} // namespace Eigen
