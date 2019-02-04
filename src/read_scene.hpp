#ifndef CCD_READ_SCENE_H
#define CCD_READ_SCENE_H

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

namespace ccd {
namespace io {
    typedef std::vector<std::vector<double>> vec2d;
    typedef std::vector<std::vector<int>> vec2i;

    void read_scene(const std::string filename, Eigen::MatrixXd& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixXd& displacements);
    void read_scene(const nlohmann::json scene, Eigen::MatrixXd& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixXd& displacements);
    void read_scene_from_str(const std::string str,  Eigen::MatrixXd& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixXd& displacements);
}
}
#endif
