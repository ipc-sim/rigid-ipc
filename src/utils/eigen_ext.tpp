#pragma once

#include "eigen_ext.hpp"

namespace Eigen {

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x)
{
    return MatrixX<T>(x.asDiagonal()).sparseView();
}

} // namespace Eigen
