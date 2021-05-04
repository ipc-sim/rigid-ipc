#pragma once

#include "serialize_json.hpp"

namespace ipc::rigid {

template <typename T, int dim, int max_dim>
void from_json(
    const nlohmann::json& json,
    Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1>& vector)
{
    typedef Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1> Vector;
    typedef std::vector<T> L;
    L list = json.template get<L>();
    vector = Eigen::Map<Vector>(list.data(), long(list.size()));
}

template <typename T>
void from_json(const nlohmann::json& json, MatrixX<T>& matrix)
{
    typedef std::vector<std::vector<T>> L;
    L list = json.get<L>();

    size_t num_rows = list.size();
    if (num_rows == 0) {
        return;
    }
    size_t num_cols = list[0].size();
    matrix.resize(long(num_rows), long(num_cols));

    for (size_t i = 0; i < num_rows; ++i) {
        assert(num_cols == list[i].size());
        matrix.row(int(i)) = Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic>>(
            list[i].data(), long(num_cols));
    }
}

template <typename T>
void from_json(const nlohmann::json& json, Matrix3<T>& matrix)
{
    typedef std::vector<std::vector<T>> L;
    L list = json.get<L>();

    size_t num_rows = list.size();
    if (num_rows == 0) {
        return;
    }
    size_t num_cols = list[0].size();
    assert(num_rows == matrix.rows() && num_cols == matrix.cols());

    for (size_t i = 0; i < num_rows; ++i) {
        assert(num_cols == list[i].size());
        matrix.row(int(i)) = Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic>>(
            list[i].data(), long(num_cols));
    }
}

template <typename Derived>
nlohmann::json to_json(const Eigen::MatrixBase<Derived>& matrix_base)
{
    using T = typename Derived::Scalar;
    size_t num_rows = size_t(matrix_base.rows());
    size_t num_cols = size_t(matrix_base.cols());

    typedef Eigen::Ref<const Eigen::Matrix<
        T, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>>
        RefMatrix;
    RefMatrix matrix(matrix_base);

    if (num_cols == 1) {
        std::vector<T> vec(matrix.data(), matrix.data() + matrix.size());
        return nlohmann::json(vec);
    }

    std::vector<std::vector<T>> mat(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        std::vector<T> row(num_cols);
        for (size_t j = 0; j < num_cols; ++j) {
            row[j] = matrix(long(i), long(j));
        }
        mat[i] = row;
    }
    return nlohmann::json(mat);
}

} // namespace ipc::rigid
