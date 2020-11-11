#pragma once

#include <Eigen/Core>
#include <nlohmann/json.hpp>

namespace swr {

template <typename T, int dim, int max_dim = dim>
void from_json(
    const nlohmann::json& json,
    Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1>& vector)
{
    typedef Eigen::Matrix<T, dim, 1, Eigen::ColMajor, max_dim, 1> Vector;
    std::vector<T> list = json.get<std::vector<T>>();
    vector = Eigen::Map<Vector>(list.data(), long(list.size()));
}

} // namespace swr
