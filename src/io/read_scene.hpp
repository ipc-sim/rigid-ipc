#ifndef CCD_READ_SCENE_H
#define CCD_READ_SCENE_H

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

namespace ccd {
namespace io {
    typedef std::vector<std::vector<double>> vec_vec_d;
    typedef std::vector<std::vector<int>> vec_vec_i;


    void read_scene(const std::string filename, Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements);
    void read_scene(const nlohmann::json scene, Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements);
    void read_scene_from_str(const std::string str,  Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements);
}
}
#endif
