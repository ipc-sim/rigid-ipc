#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

namespace ccd {
namespace io {

    template <typename T>
    nlohmann::json to_json(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);

    template <typename T>
    void from_json(const nlohmann::json&,
        Eigen::Matrix<T, Eigen::Dynamic, 1>& vector);

    template <typename T>
    nlohmann::json to_json(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix);

    template <typename type>
    void from_json(const nlohmann::json&,
        Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>& matrix);

    template <>
    void from_json<bool>(const nlohmann::json& json,
        Eigen::Matrix<bool, Eigen::Dynamic, 1>& vector);

} // namespace io
} // namespace ccd

#include "serialize_json.tpp"
