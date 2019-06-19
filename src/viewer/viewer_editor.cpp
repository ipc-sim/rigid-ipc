#include "viewer.hpp"

namespace ccd {

// Select all vertices (clears selected displacements).
void ViewerMenu::select_all_vertices()
{
    state.selected_displacements.clear();
    Eigen::VectorXi all_idxs = Eigen::VectorXi::LinSpaced(
        state.vertices.rows(), 0, state.vertices.rows() - 1);
    state.selected_points
        = std::vector<int>(all_idxs.data(), all_idxs.data() + all_idxs.size());
}

// Select all displacments (clears selected vertices).
void ViewerMenu::select_all_displacements()
{
    state.selected_points.clear();
    Eigen::VectorXi all_idxs = Eigen::VectorXi::LinSpaced(
        state.displacements.rows(), 0, state.displacements.rows() - 1);
    state.selected_displacements
        = std::vector<int>(all_idxs.data(), all_idxs.data() + all_idxs.size());
}

// Select the vertices of displacements that are connected to the selection.
void ViewerMenu::select_connected()
{
    // Determine which selection we want to extend
    std::vector<int>& selection = state.selected_points.size() > 0
        ? state.selected_points
        : state.selected_displacements;
    // Build an adjacency list representation of the mesh
    const auto adjacency_list = state.create_adjacency_list();

    std::unordered_set<int> new_selection;
    for (const auto& selected_id : selection) {
        state.find_connected_vertices(
            selected_id, adjacency_list, new_selection);
    }
    selection.assign(new_selection.begin(), new_selection.end());

    recolor_edges();
    recolor_displacements();
}

// Duplicate selected vertices and edges that have both end-points selected.
void ViewerMenu::duplicate_selected()
{
    long num_old_vertices = state.vertices.rows();
    Eigen::Vector2d delta_com(0, -1);
    state.duplicate_selected_vertices(delta_com);
    state.selected_points.clear();

    state_history.push_back(state);
    load_state();

    // Select the newly created points
    for (long i = num_old_vertices; i < state.vertices.rows(); i++) {
        state.selected_points.push_back(i);
    }
    recolor_edges();
}

// Loop over edges creating a new vertex in the center.
void ViewerMenu::subdivide_edges()
{
    // Save the number of old and new vertices
    const long num_old_vertices = state.vertices.rows();
    const long num_new_vertices = state.edges.rows();

    // Resize the fields
    state.vertices.conservativeResize(
        num_old_vertices + num_new_vertices, state.vertices.cols());
    state.displacements.conservativeResize(
        num_old_vertices + num_new_vertices, state.displacements.cols());
    state.edges.conservativeResize(2 * state.edges.rows(), 2);
    // Copy the fixed dof because it is reset during reset_scene
    auto is_dof_fixed = state.getOptimizationProblem().is_dof_fixed;
    is_dof_fixed.resize(num_old_vertices, 2);
    is_dof_fixed.conservativeResize(state.displacements.rows(), 2);

    for (long i = 0; i < num_new_vertices; i++) {
        // Average adjacent vertices.
        state.vertices.row(num_old_vertices + i) = 0.5
            * (state.vertices.row(state.edges(i, 1))
                + state.vertices.row(state.edges(i, 0)));

        // Average adjacent displacements.
        state.displacements.row(num_old_vertices + i) = 0.5
            * (state.displacements.row(state.edges(i, 1))
                + state.displacements.row(state.edges(i, 0)));

        // AND endpoints of edge to determine if the new vertex is fixed.
        for (int j = 0; j < 2; j++) {
            is_dof_fixed(num_old_vertices + i, j)
                = is_dof_fixed(state.edges(i, 0), j)
                && is_dof_fixed(state.edges(i, 1), j);
        }

        // Split the edge into two edges.
        state.edges.row(num_new_vertices + i) = state.edges.row(i);
        state.edges(i, 1) = num_old_vertices + i;
        state.edges(num_new_vertices + i, 0) = num_old_vertices + i;
    }

    state.reset_scene();

    // Flatten the is_dof_fixed
    is_dof_fixed.resize(is_dof_fixed.size(), 1);
    state.getOptimizationProblem().is_dof_fixed = is_dof_fixed;

    state_history.push_back(state);
    load_state();
}

// Smooth the vertices using a weighted average of the adjacent vertices.
void ViewerMenu::smooth_vertices()
{
    // Weight for the original position
    double w0 = 0.75; // TODO: Expose this parameter

    // Copy the vertices/displacements over and weight them
    Eigen::MatrixX2d smoothed_vertices = w0 * state.vertices;
    Eigen::MatrixX2d smoothed_displacements = w0 * state.displacements;

    auto adjacency_list = state.create_adjacency_list();

    for (long i = 0; i < state.vertices.rows(); i++) {
        const unsigned long num_neighbors = adjacency_list[i].size();
        for (const auto& vertex_idx : adjacency_list[i]) {
            double w_adjacent = (1 - w0) / num_neighbors;
            // Equally weight adjacent vertex positions
            smoothed_vertices.row(i)
                += w_adjacent * state.vertices.row(vertex_idx);
            // Equally weight adjacent displacements
            smoothed_displacements.row(i)
                += w_adjacent * state.displacements.row(vertex_idx);
        }
    }

    // Update the state
    state.vertices = smoothed_vertices;
    state.displacements = smoothed_displacements;

    auto is_dof_fixed = state.getOptimizationProblem()
                            .is_dof_fixed; // Save this from being reset
    state.reset_scene();
    state.getOptimizationProblem().is_dof_fixed = is_dof_fixed;

    state_history.push_back(state);
    load_state();
}

} // namespace ccd
