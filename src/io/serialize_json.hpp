#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

template <typename T, int dim, int max_dim = dim>
void from_json(
    const nlohmann::json& json,
    Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1>& vector);

template <typename T>
void from_json(const nlohmann::json& json, MatrixX<T>& matrix);

template <typename T>
void from_json(const nlohmann::json& json, Matrix3<T>& matrix);

template <>
void from_json<bool>(const nlohmann::json& json, VectorX<bool>& vector);

template <typename Derived>
nlohmann::json to_json(const Eigen::MatrixBase<Derived>& matrix);

nlohmann::json to_json_string(
    const Eigen::MatrixXd& matrix, const std::string& format = ".16e");

} // namespace ipc::rigid

#include "serialize_json.tpp"
