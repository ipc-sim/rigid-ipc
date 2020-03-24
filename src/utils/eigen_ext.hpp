#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Eigen {

typedef Array<bool, Dynamic, 1> ArrayXb;
typedef Matrix<bool, Dynamic, 1> VectorXb;
typedef Matrix<bool, 3, 1> Vector3b;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

template <typename T, int dim, int max_dim = dim>
using Vector = Matrix<T, dim, 1, ColMajor, max_dim, 1>;
template <typename T> using Vector1 = Vector<T, 1>;
template <typename T> using Vector2 = Vector<T, 2>;
template <typename T> using Vector3 = Vector<T, 3>;
template <typename T> using VectorX = Vector<T, Eigen::Dynamic>;
template <typename T> using Matrix2 = Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Matrix<T, 3, 3>;
template <typename T> using MatrixX = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
template <typename T> using VectorX3 = Vector<T, Eigen::Dynamic, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
using VectorX3d = VectorX3<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 6 x 1
template <typename T> using VectorX6 = Vector<T, Eigen::Dynamic, 6>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
using VectorX6d = VectorX6<double>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 1
using VectorX6b = VectorX6<bool>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
template <typename T>
using MatrixXX3 = Matrix<T, Eigen::Dynamic, Eigen::Dynamic, ColMajor, 3, 3>;
/// @brief A dynamic size matrix with a fixed maximum size of 3 x 3
using MatrixXX3d = MatrixXX3<double>;

template <typename T> SparseMatrix<T> SparseDiagonal(const VectorX<T>& x);

template <typename T> MatrixXX3<T> Hat(Eigen::VectorX3<T> x);

} // namespace Eigen

#include "eigen_ext.tpp"
