#include "eigen_ext.hpp"

namespace Eigen {

template <typename T>
SparseMatrix<T> SparseDiagonal(const Matrix<T, Dynamic, 1>& x)
{
    return Matrix<T, Dynamic, Dynamic>(x.asDiagonal()).sparseView();
}

template SparseMatrix<double> SparseDiagonal(
    const Matrix<double, Dynamic, 1>& x);

} // namespace Eigen
