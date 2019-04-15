#include "viewer.hpp"

#include <igl/opengl/ViewerData.h>

namespace ccd {

// ---------------------------------------------------------------------
// HELPER FUNCTIONS
// ---------------------------------------------------------------------
void viewer_set_vertices(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& V);
void viewer_set_points(igl::opengl::ViewerData& data, const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& C);
void viewer_set_edges(igl::opengl::ViewerData& data, const Eigen::MatrixXd& P,
    const Eigen::MatrixXi& E, const Eigen::MatrixXd& C);
void viewer_set_edges2(igl::opengl::ViewerData& data, const Eigen::MatrixXd& P1,
    const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C);
void viewer_update_edges(igl::opengl::ViewerData& data,
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges);
void viewer_update_edges2(igl::opengl::ViewerData& data,
    const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2);
void viewer_update_points(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& vertices);
void viewer_update_vertices(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& vertices);
void viewer_color_edges(
    igl::opengl::ViewerData& data, const Eigen::RowVector3d color);
void viewer_color_points(
    igl::opengl::ViewerData& data, const Eigen::RowVector3d color);
void viewer_highlight_points(igl::opengl::ViewerData& data,
    const std::vector<int>& nodes, const Eigen::RowVector3d& color_hl);
void viewer_highlight_edge(igl::opengl::ViewerData& data, const int edge,
    const Eigen::RowVector3d& color_hl);

// ---------------------------------------------------------------------
// VIEWER DRAW FUNCTIONS
// ---------------------------------------------------------------------

void ViewerMenu::redraw_scene()
{
    redraw_at_time();
    redraw_displacements();
    redraw_opt_displacements();
    redraw_grad_volume(state.use_opt_gradient);
}

void ViewerMenu::create_edges()
{
    auto& data = viewer->data_list[edges_data_id];
    viewer_set_points(data, state.vertices, color_edge);
    viewer_set_edges(data, state.vertices, state.edges, color_edge);
    viewer_set_vertices(data, state.vertices);
    data.point_size = 10 * pixel_ratio();
}

void ViewerMenu::redraw_at_time() {
    redraw_edges(state.get_vertex_at_time());
    redraw_grad_volume(/*use_opt_volume=*/false);
}

void ViewerMenu::redraw_at_opt_time()
{
    redraw_edges(state.get_opt_vertex_at_time());
    redraw_grad_volume(/*use_opt_volume=*/true);
}

void ViewerMenu::redraw_edges(const Eigen::MatrixXd& vertices)
{
    auto& data = viewer->data_list[edges_data_id];
    viewer_update_edges(data, vertices, state.edges);
    viewer_update_points(data, vertices);
    viewer_update_vertices(data, vertices);
}

void ViewerMenu::recolor_edges()
{
    auto& data = viewer->data_list[edges_data_id];
    viewer_color_edges(data, color_edge);
    viewer_color_points(data, color_edge);

    viewer_highlight_points(data, state.selected_points, color_sl);
}

void ViewerMenu::create_displacements()
{
    auto& data = viewer->data_list[displ_data_id];
    const Eigen::MatrixXd& P2 = state.vertices + state.displacements;
    viewer_set_points(data, P2, color_displ);
    viewer_set_edges2(data, state.vertices, P2, color_displ);
    data.point_size = 1 * pixel_ratio();
}

void ViewerMenu::redraw_displacements()
{
    const Eigen::MatrixXd& P2 = state.vertices + state.displacements;
    viewer_update_edges2(
        viewer->data_list[displ_data_id], state.vertices, P2);
    viewer_update_points(viewer->data_list[displ_data_id], state.vertices + state.displacements);
}

void ViewerMenu::recolor_displacements()
{
    auto& data = viewer->data_list[displ_data_id];
    viewer_color_edges(data, color_displ);
    viewer_color_points(data, color_displ);
    viewer_highlight_points(data, state.selected_displacements, color_sl);
}

void ViewerMenu::create_opt_displacements()
{
    auto& data = viewer->data_list[opt_displ_data_id];
    viewer_set_edges2(data, state.vertices, state.vertices, color_opt_displ);
}

void ViewerMenu::redraw_opt_displacements()
{
    const Eigen::MatrixXd& P2 = state.vertices + state.get_opt_displacements();
    viewer_update_edges2(
        viewer->data_list[opt_displ_data_id], state.vertices, P2);
}

void ViewerMenu::recolor_opt_displacements()
{
    auto& data = viewer->data_list[opt_displ_data_id];
    viewer_color_edges(data, color_opt_displ);
}

void ViewerMenu::create_grad_volume()
{
    auto& data = viewer->data_list[gradient_data_id];
    viewer_set_edges2(data, state.vertices, state.vertices, color_grad);
}

void ViewerMenu::redraw_grad_volume(const bool use_opt_gradient)
{
    state.use_opt_gradient = use_opt_gradient;
    Eigen::MatrixX2d V, grad;

    if (use_opt_gradient) {
        V = state.get_opt_vertex_at_time();
        grad = state.get_opt_volume_grad().normalized();

    } else {
        V = state.get_vertex_at_time();
        grad = state.get_volume_grad().normalized();
    }

    auto& data = viewer->data_list[gradient_data_id];
    viewer_update_edges2(data, V, V + grad * double(state.grad_scaling));
}

void ViewerMenu::recolor_grad_volume()
{
    viewer_color_edges(viewer->data_list[gradient_data_id], color_grad);
}
// ---------------------------------------------------------------------
// HELPER FUNCTIONS IMPLEMENTATION
// ---------------------------------------------------------------------
void viewer_set_vertices(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& V)
{
    Eigen::MatrixXd V_temp;

    // If P only has two columns, pad with a column of zeros
    if (V.cols() == 2) {
        V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
        V_temp.block(0, 0, V.rows(), 2) = V;
    } else {
        V_temp = V;
    }
    data.set_vertices(V_temp);

    Eigen::MatrixXd nz = Eigen::MatrixXd::Zero(V.rows(), 3);
    nz.col(2).setConstant(1.0);
    data.set_normals(nz);
}

void viewer_set_points(igl::opengl::ViewerData& data, const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& C)
{
    Eigen::MatrixXd V_temp;

    // If P only has two columns, pad with a column of zeros
    if (V.cols() == 2) {
        V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
        V_temp.block(0, 0, V.rows(), 2) = V;
    } else {
        V_temp = V;
    }
    data.set_points(V_temp, C);
}

void viewer_update_points(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& vertices)
{
    long dim = vertices.cols();

    Eigen::MatrixXd& points = data.points;
    for (unsigned i = 0; i < points.rows(); ++i) {
        points.row(i).segment(0, dim) = vertices.row(i);
    }

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
}

void viewer_update_vertices(
    igl::opengl::ViewerData& data, const Eigen::MatrixXd& V)
{
    Eigen::MatrixXd V_temp;

    // If P only has two columns, pad with a column of zeros
    if (V.cols() == 2) {
        V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
        V_temp.block(0, 0, V.rows(), 2) = V;
    } else {
        V_temp = V;
    }
    data.set_vertices(V_temp);
}

void viewer_set_edges(igl::opengl::ViewerData& data, const Eigen::MatrixXd& P,
    const Eigen::MatrixXi& E, const Eigen::MatrixXd& C)
{
    Eigen::MatrixXd P_temp;

    // If P only has two columns, pad with a column of zeros
    if (P.cols() == 2) {
        P_temp = Eigen::MatrixXd::Zero(P.rows(), 3);
        P_temp.block(0, 0, P.rows(), 2) = P;
    } else {
        P_temp = P;
    }
    data.set_edges(P_temp, E, C);
}

void viewer_set_edges2(igl::opengl::ViewerData& data, const Eigen::MatrixXd& P1,
    const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C)
{
    data.lines.resize(0, 9);
    data.add_edges(P1, P2, C);
}

void viewer_update_edges(igl::opengl::ViewerData& data,
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
{
    long dim = vertices.cols();

    Eigen::MatrixXd& lines = data.lines;
    assert(data.lines.rows() == edges.rows());

    for (unsigned i = 0; i < edges.rows(); ++i) {
        lines.row(i).segment(0, dim) << vertices.row(edges(i, 0));
        lines.row(i).segment(3, dim) << vertices.row(edges(i, 1));
    }

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
}

void viewer_update_edges2(igl::opengl::ViewerData& data,
    const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2)
{
    long dim = P1.cols();

    Eigen::MatrixXd& lines = data.lines;
    assert(data.lines.rows() == P1.rows());

    for (unsigned i = 0; i < P1.rows(); ++i) {
        lines.row(i).segment(0, dim) << P1.row(i);
        lines.row(i).segment(3, dim) << P2.row(i);
    }

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
}

void viewer_color_edges(
    igl::opengl::ViewerData& data, const Eigen::RowVector3d color)
{
    Eigen::MatrixXd& lines = data.lines;
    lines.block(0, 6, lines.rows(), 3).rowwise() = color;

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
}

void viewer_color_points(
    igl::opengl::ViewerData& data, const Eigen::RowVector3d color)
{
    Eigen::MatrixXd& points = data.points;
    points.block(0, 3, points.rows(), 3).rowwise() = color;
    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
}

void viewer_highlight_points(igl::opengl::ViewerData& data,
    const std::vector<int>& nodes, const Eigen::RowVector3d& color_hl)
{
    Eigen::MatrixXd& points = data.points;
    for (unsigned i = 0; i < nodes.size(); ++i) {
        points.row(nodes[i]).segment(3, 3) = color_hl;
    }

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
}

void viewer_highlight_edge(igl::opengl::ViewerData& data, const int edge,
    const Eigen::RowVector3d& color_hl)
{
    Eigen::MatrixXd& lines = data.lines;
    lines.row(edge).segment(6, 3) = color_hl;

    data.dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
}

} // namespace ccd
