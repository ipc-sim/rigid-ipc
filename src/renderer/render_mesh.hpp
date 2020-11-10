#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include <nlohmann/json.hpp>

#include <renderer/attributes.hpp>

struct Scene {
    Eigen::Vector4F background_color;
    Eigen::Vector3F ambient_light;
    Float line_thickness;

    Camera camera;
    std::vector<Material> materials;
    std::vector<Light> lights;

    Scene(const nlohmann::json&);
};

/// Render a mesh with vertices V, edges E, faces F, and vertex colors C
bool render_mesh(
    const Scene& scene,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& C,
    const std::string& filename);
