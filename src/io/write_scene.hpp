#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

namespace ccd {
namespace io {
    template <typename type>
    std::vector<std::vector<type>> to_vec2(
        const Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>& matrix);

    void write_scene(const std::string filename,
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges,
        const Eigen::MatrixXd& displacements);
    nlohmann::json write_scene(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges, const Eigen::MatrixXd& displacements);
} // namespace io
} // namespace ccd

