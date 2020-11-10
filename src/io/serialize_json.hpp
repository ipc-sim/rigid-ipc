#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace io {

    template <typename T>
    nlohmann::json to_json(const Eigen::VectorX<T>& vector);

    template <typename T, int dim, int max_dim = dim>
    void from_json(
        const nlohmann::json& json,
        Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1>& vector);

    template <typename T>
    nlohmann::json to_json(const Eigen::MatrixX<T>& matrix);

    nlohmann::json to_json_string(
        const Eigen::MatrixXd& matrix, const std::string& format = ".16e");

    template <typename T>
    void from_json(const nlohmann::json& json, Eigen::MatrixX<T>& matrix);

    template <>
    void
    from_json<bool>(const nlohmann::json& json, Eigen::VectorX<bool>& vector);

} // namespace io
} // namespace ccd

#include "serialize_json.tpp"
