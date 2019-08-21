#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Eigen {
typedef Array<bool, Dynamic, 1> ArrayXb;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<bool, 3, 1> Vector3b;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

template <typename T>
SparseMatrix<T> SparseDiagonal(const Matrix<T, Dynamic, 1>& x);

} // namespace Eigen

#include "eigen_ext.tpp"
