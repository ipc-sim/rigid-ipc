#include "viewer.hpp"
#include <igl/slice.h>
#include <igl/unproject_onto_mesh.h>

namespace ccd {

// ----------------------------------------------------------------------------
// HELPER FUNCTIONS
// ----------------------------------------------------------------------------

bool update_selection(const Eigen::MatrixXd& points, const int modifier,
    const Eigen::RowVector2d& click_coord, const double thr,
    std::vector<int>& selection)
{
    // If a vertex was clicked while holding "shift", it is added to the
    // selected points. If "shift" was not clicked, the selection is replaced by
    // the clicked-vertex.
    int index;
    double dist
        = (points.rowwise() - click_coord).rowwise().norm().minCoeff(&index);
    if (dist < thr) {
        auto index_position
            = std::find(selection.begin(), selection.end(), index);
        if (index_position != selection.end()) {
            selection.erase(index_position);
            return false;
        } else if (modifier == GLFW_MOD_SHIFT) {
            selection.push_back(index);
        } else {
            selection.clear();
            selection.push_back(index);
        }
    } else {
        selection.clear();
        return false;
    }
    return true;
}

void translation_delta(const Eigen::MatrixXd& points, const int modifier,
    const Eigen::RowVector2d& click_coord, const std::vector<int>& selection,
    Eigen::RowVector2d& delta)
{
    // Compute the center of mass of the selecteion
    Eigen::MatrixXd selected_points(selection.size(), points.cols());
    igl::slice(points, Eigen::VectorXi::Map(selection.data(), selection.size()),
        Eigen::VectorXi::LinSpaced(
            selected_points.cols(), 0, selected_points.cols() - 1),
        selected_points);
    Eigen::RowVector2d center
        = selected_points.colwise().sum() / selected_points.rows();

    // Translate the center of mass to the clicked point
    delta = click_coord - center;

    // If shift is pressed, zero out all dimensions except for the maxium one.
    if (modifier == GLFW_MOD_SHIFT) {
        int index;
        delta.cwiseAbs().maxCoeff(&index); // Get the index of the max coeff

        double aux = delta[index]; // Save max coeff
        delta.setZero();           // Zero out delta
        delta[index] = aux;        // Restore max coeff
    }
}

// ----------------------------------------------------------------------------
// CLASS FUNCTIONS
// ----------------------------------------------------------------------------

void ViewerMenu::undo()
{
    if (state_history.size() > 1) {
        state_history.pop_back();
        state = state_history.back();
        load_state();
    }
}

void ViewerMenu::clicked__select(
    const int /*button*/, const int modifier, const Eigen::RowVector2d& coord)
{
    double thr = state.canvas_width / 100;
    bool clicked_vertex = update_selection(
        state.vertices, modifier, coord, thr, state.selected_points);
    if (clicked_vertex) {
        state.selected_displacements.clear();
    } else {
        update_selection(state.vertices + state.displacements, modifier, coord,
            thr, state.selected_displacements);
    }
    recolor_edges();
    recolor_displacements();
}

void ViewerMenu::clicked__translate(
    const int /*button*/, const int modifier, const Eigen::RowVector2d& coord)
{
    if (state.selected_points.size() > 0) {
        Eigen::RowVector2d delta;
        translation_delta(
            state.vertices, modifier, coord, state.selected_points, delta);
        for (int selected_idx : state.selected_points) {
            state.move_vertex(selected_idx, delta);
        }
    } else if (state.selected_displacements.size() > 0) {
        Eigen::RowVector2d delta;
        translation_delta(state.vertices + state.displacements, modifier, coord,
            state.selected_displacements, delta);
        for (int selected_idx : state.selected_displacements) {
            state.move_displacement(selected_idx, delta);
        }
    }

    state_history.push_back(state);
    redraw_scene();
}
void ViewerMenu::clicked__add_node(const int /*button*/, const int /*modifier*/,
    const Eigen::RowVector2d& coord)
{
    state.add_vertex(coord);
    state_history.push_back(state);
    load_state();
}

void ViewerMenu::connect_selected_vertices()
{
    size_t num_sl = state.selected_points.size();
    if (num_sl < 2) {
        return;
    }
    Eigen::MatrixXi new_edges(num_sl - 1, 2);
    for (uint i = 0; i < num_sl - 1; ++i) {
        new_edges.row(i) << state.selected_points[i],
            state.selected_points[i + 1];
    }
    state.add_edges(new_edges);
    state_history.push_back(state);
    load_state();
}

void ViewerMenu::clicked_on_canvas(
    const int button, const int modifier, const Eigen::RowVector2d& coord)
{
    switch (edit_mode) {
    case select:
        clicked__select(button, modifier, coord);
        break;
    case translate:
        clicked__translate(button, modifier, coord);
        break;
    case add_node:
        clicked__add_node(button, modifier, coord);
        break;
    default:
        break;
    }
}

bool ViewerMenu::mouse_down(int button, int modifier)
{
    if (Super::mouse_down(button, modifier)) {
        return true;
    }

    // pick a node
    double x = viewer->current_mouse_x;
    double y = double(viewer->core().viewport(3)) - viewer->current_mouse_y;

    int face_id;
    Eigen::Vector3d barycentric_coords;

    Eigen::MatrixXd BV = viewer->data_list[canvas_data_id].V;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer->core().view,
            viewer->core().proj, viewer->core().viewport, BV, canvas_faces,
            face_id, barycentric_coords)) {
        Eigen::Vector3i face_nodes = canvas_faces.row(face_id);
        Eigen::RowVector3d coords
            = BV.row(face_nodes[0]) * barycentric_coords[0]
            + BV.row(face_nodes[1]) * barycentric_coords[1]
            + BV.row(face_nodes[2]) * barycentric_coords[2];

        // pass 2d coords
        clicked_on_canvas(button, modifier, coords.segment(0, 2));
    }
    return true;
}

bool ViewerMenu::key_pressed(unsigned int key, int modifiers)
{
    if (Super::key_pressed(key, modifiers)) {
        return true;
    }

    if (key == 'q') {
        edit_mode = ViewerEditMode::select;
    } else if (key == 'w') {
        edit_mode = ViewerEditMode::translate;
    } else if (key == 'Z' && modifiers == GLFW_MOD_SHIFT) {
        undo();
    } else {
        return false;
    }
    return true;
}

} // namespace ccd
