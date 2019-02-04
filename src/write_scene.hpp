#ifndef CCD_WRITE_SCENE_HPP
#define CCD_WRITE_SCENE_HPP

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

namespace ccd {
namespace io {
    void write_scene(const std::string filename, const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges, const Eigen::MatrixXd& displacements);
    nlohmann::json write_scene(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges, const Eigen::MatrixXd& displacements);
}
}
#endif
