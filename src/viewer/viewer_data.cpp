#include "viewer_data.hpp"

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include <igl/remove_unreferenced.h>
#include <igl/slice.h>

#include <iostream>

namespace ipc::rigid {

// =============================================================================
/// MeshData
// =============================================================================

const glm::vec3 MeshData::EDGE_COLOR(1.0, 1.0, 1.0); // #FFFFFF
const glm::vec3
    MeshData::STATIC_COLOR(0xB3 / 255.0, 0xB3 / 255.0, 0xB3 / 255.0); // #B3B3B3
// const glm::vec3 MeshData::KINEMATIC_COLOR(1.0, 0.5, 0.0); // #FF8000
const glm::vec3 MeshData::KINEMATIC_COLOR(
    0xE3 / 255.0, 0x82 / 255.0, 0x1C / 255.0); // #E3821C

void MeshData::set_mesh(
    const Eigen::MatrixXd& _V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<int>& CE2E,
    const std::vector<int>& CV2V)
{
    const int dim = _V.cols();
    Eigen::MatrixXd V(_V.rows(), 3);
    V.leftCols(dim) = _V;
    if (dim == 2) {
        V.col(2).setZero();
    }

    if (F.size() > 0) {
        ps_surface_mesh = ps::registerSurfaceMesh("mesh", V, F);
        ps_velocities = ps_surface_mesh->addVertexVectorQuantity(
            "velocities", Eigen::MatrixXd::Zero(V.rows(), V.cols()));
        ps_velocities->setEnabled(false);
        ps_velocities->setVectorColor(velocity_color);
    }

    if (CE2E.size() > 0) {
        Eigen::MatrixXi CE;
        igl::slice(
            E, Eigen::Map<const Eigen::VectorXi>(CE2E.data(), CE2E.size()), 1,
            CE);
        assert(CE.rows() == CE2E.size() && CE.cols() == E.cols());

        Eigen::MatrixXd CE_V;
        Eigen::MatrixXi CE_E;
        Eigen::VectorXi CE_I;
        igl::remove_unreferenced(V, CE, CE_V, CE_E, CE_I, CE_J);

        if (CE.rows() != 0) {
            ps_codim_edges = ps::registerCurveNetwork("edges", CE_V, CE_E);
            ps_codim_edges->setRadius(dim == 2 ? 1e-3 : 5e-4);
            if (dim == 3) {
                ps_codim_edges->setColor(EDGE_COLOR);
            }
        }
    } else if (ps_codim_edges != nullptr) {
        ps::removeStructure(ps_codim_edges);
        ps_codim_edges = nullptr;
    }

    if (CV2V.size() > 0) {
        CV_J = Eigen::Map<const Eigen::VectorXi>(CV2V.data(), CV2V.size());

        Eigen::MatrixXd CV;
        igl::slice(V, CV_J, 1, CV);

        ps_codim_points = ps::registerPointCloud("vertices", CV);
        ps_codim_points->setPointRadius(5e-4);
        ps_codim_points->setPointColor(EDGE_COLOR);
    } else if (ps_codim_points != nullptr) {
        ps::removeStructure(ps_codim_points);
        ps_codim_points = nullptr;
    }
}

void MeshData::update_vertices(const Eigen::MatrixXd& _V)
{
    const int dim = _V.cols();
    Eigen::MatrixXd V(_V.rows(), 3);
    V.leftCols(dim) = _V;
    if (dim == 2) {
        V.col(2).setZero();
    }

    if (ps_surface_mesh != nullptr) {
        ps_surface_mesh->updateVertexPositions(V);
    }
    if (ps_codim_edges != nullptr) {
        Eigen::MatrixXd CE_V;
        igl::slice(V, CE_J, 1, CE_V);
        ps_codim_edges->updateNodePositions(CE_V);
    }
    if (ps_codim_points != nullptr) {
        Eigen::MatrixXd CV;
        igl::slice(V, CV_J, 1, CV);
        ps_codim_points->updatePointPositions(CV);
    }
}

void MeshData::update_velocities(const Eigen::MatrixXd& velocities)
{
    if (ps_surface_mesh != nullptr) {
        ps_velocities =
            ps_surface_mesh->addVertexVectorQuantity("velocities", velocities);
    }
}

void MeshData::set_vertex_types(const Eigen::VectorXi& vertex_types)
{
    std::vector<glm::vec3> colors(vertex_types.size());
    for (int i = 0; i < vertex_types.size(); i++) {
        if (vertex_types(i) == 0) {
            colors[i] = STATIC_COLOR;
        } else if (vertex_types(i) == 1) {
            colors[i] = KINEMATIC_COLOR;
        } else {
            colors[i] = mesh_color; // Default color for other types
        }
    }
    if (ps_surface_mesh != nullptr) {
        ps_surface_mesh->addVertexColorQuantity("type", colors);
    } else if (ps_codim_edges != nullptr) {
        ps_codim_edges->addNodeColorQuantity("type", colors);
    }
}

void MeshData::set_vertex_ids(const Eigen::VectorXi& ids)
{
    if (ps_surface_mesh != nullptr) {
        auto tmp = ps_surface_mesh->addVertexScalarQuantity("id", ids);
        tmp->setEnabled(true);
        tmp->setColorMap("rainbow");
    } else if (ps_codim_edges != nullptr) {
        auto tmp = ps_codim_edges->addNodeScalarQuantity("id", ids);
        tmp->setEnabled(true);
        tmp->setColorMap("rainbow");
    }
}

// =============================================================================
/// CoMData
// =============================================================================

void CoMData::set_coms(const ipc::rigid::PosesD& poses)
{
    int dim = poses.size() ? poses[0].dim() : 0;
    Eigen::MatrixXd coms = Eigen::MatrixXd::Zero(poses.size(), 3);
    std::vector<Eigen::MatrixXd> principle_axes(dim, coms);

    for (int i = 0; i < poses.size(); i++) {
        coms.row(i).head(dim) = poses[i].position;
        MatrixMax3d R = poses[i].construct_rotation_matrix();
        for (int j = 0; j < dim; j++) {
            VectorMax3d axis = VectorMax3d::Zero(dim);
            axis(j) = 1;
            principle_axes[j].row(i).head(dim) = R * axis;
        }
    }

    ps_coms = ps::registerPointCloud("CoMs", coms);
    ps_coms->setPointColor(com_color);

    for (int j = 0; j < dim; j++) {
        auto ps_axis = ps_coms->addVectorQuantity(
            fmt::format("{}", "xyz"[j]), principle_axes[j]);
        ps_axis->setEnabled(true);
        ps_axis->setVectorColor(
            j == 0
                ? glm::vec3(1.0, 0.14902, 0.0)                      // #FF2600
                : (j == 1 ? glm::vec3(0.0, 0.97647, 0.0)            // #00F900
                          : glm::vec3(0.01176, 0.16078, 0.81569))); // #0329D0
    }
}

void CoMData::update_coms(const ipc::rigid::PosesD& poses)
{
    int dim = poses.size() ? poses[0].dim() : 0;
    Eigen::MatrixXd coms = Eigen::MatrixXd::Zero(poses.size(), 3);
    std::vector<Eigen::MatrixXd> principle_axes(dim, coms);

    for (int i = 0; i < poses.size(); i++) {
        coms.row(i).head(dim) = poses[i].position;
        MatrixMax3d R = poses[i].construct_rotation_matrix();
        for (int j = 0; j < dim; j++) {
            VectorMax3d axis = VectorMax3d::Zero(dim);
            axis(j) = 1;
            principle_axes[j].row(i).head(dim) = R * axis;
        }
    }

    ps_coms = ps::registerPointCloud("CoMs", coms);

    for (int j = 0; j < dim; j++) {
        ps_coms->addVectorQuantity(
            fmt::format("{}", "xyz"[j]), principle_axes[j]);
    }
}

} // namespace ipc::rigid
