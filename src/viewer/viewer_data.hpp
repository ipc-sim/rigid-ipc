#pragma once

#include <physics/pose.hpp>

#include <Eigen/Core>

#include <glm/glm.hpp>

namespace polyscope {
class SurfaceMesh;
class PointCloud;
class CurveNetwork;
class SurfaceVertexVectorQuantity;
} // namespace polyscope
namespace ps = polyscope;

namespace ipc::rigid {

class MeshData {
public:
    MeshData() = default;

    void set_mesh(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const std::vector<int>& CE2E,
        const std::vector<int>& CV2V);

    void update_vertices(const Eigen::MatrixXd& V);
    void update_velocities(const Eigen::MatrixXd& velocities);
    void set_vertex_types(const Eigen::VectorXi& vertex_types);
    void set_vertex_ids(const Eigen::VectorXi& ids);

    glm::vec3 mesh_color = glm::vec3(
        0xE3 / 255.0, 0x1C / 255.0, 0x1C / 255.0); // #E31C1C - ALIZARIN
    // 0xE7 / 255.0, 0x4C / 255.0, 0x3C / 255.0); // #E74C3C - ALIZARIN

    glm::vec3 velocity_color = glm::vec3(
        0xF1 / 255.0, 0xC4 / 255.0, 0x0F / 255.0); // #F1C40F - SUN FLOWER

private:
    ps::SurfaceMesh* ps_surface_mesh = nullptr;
    ps::CurveNetwork* ps_codim_edges = nullptr;
    ps::PointCloud* ps_codim_points = nullptr;
    ps::SurfaceVertexVectorQuantity* ps_velocities = nullptr;

    Eigen::VectorXi m_vertex_type;

    Eigen::VectorXi CE_J;
    Eigen::VectorXi CV_J;

    static const glm::vec3 EDGE_COLOR;
    static const glm::vec3 STATIC_COLOR;
    static const glm::vec3 KINEMATIC_COLOR;

    std::vector<std::string> vertex_data_labels;
};

class CoMData {
public:
    CoMData() = default;

    void set_coms(const ipc::rigid::PosesD& poses);
    void update_coms(const ipc::rigid::PosesD& poses);

    glm::vec3 com_color = glm::vec3(
        0xF1 / 255.0, 0xC4 / 255.0, 0x0F / 255.0); // #F1C40F - SUN FLOWER

private:
    ps::PointCloud* ps_coms = nullptr;
};

} // namespace ipc::rigid
