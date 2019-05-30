// Draw the viewer menu buttons and fields.

#include "viewer.hpp"

#include <unordered_set>

#include <logger.hpp>
#include <viewer/constraint_view.hpp>
#include <viewer/imgui_ext.hpp>
#include <viewer/solver_view.hpp>

namespace ccd {

/**
 * @brief double_color_edit: util function to draw a colorwheel
 * @param label
 * @param color
 */
bool double_color_edit(const char* label, Eigen::RowVector3d& color)
{
    Eigen::Vector3f color_f = color.cast<float>();
    bool changed = false;
    if (ImGui::ColorEdit3(label, color_f.data(),
            ImGuiColorEditFlags_NoInputs
                | ImGuiColorEditFlags_PickerHueWheel)) {
        color = color_f.cast<double>();
        changed = true;
    }
    return changed;
}

void color_button()
{
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5, 0.0, 0.0, 1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f, 0.0, 0.0, 1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(1.0, 0.0, 0.0, 1.0));
}

void ViewerMenu::draw_menu()
{
    float menu_width = 200.f * menu_scaling();
    draw_labels_window();

    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _viewer_menu_visible = true;

    // ---------------------------------------------------------------------------------
    ImGui::Begin("CCD Viewer", &_viewer_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);

    draw_io();

    const char* level_strings[] = SPDLOG_LEVEL_NAMES;
    if (ImGui::Combo("log-level##logger", &state.log_level, level_strings,
            CCD_IM_ARRAYSIZE(level_strings))) {
        spdlog::set_level(
            static_cast<spdlog::level::level_enum>(state.log_level));
    }

    ImGui::Separator();
    ImGui::Separator();
    draw_ccd_steps();

    ImGui::Separator();
    ImGui::Separator();
    draw_optimization();
    draw_optimization_results();

    ImGui::Separator();
    ImGui::Separator();
    draw_edit_modes();
    ImGui::PopItemWidth();
    ImGui::End();

    // ---------------------------------------------------------------------------------
    ImGui::SetNextWindowPos(
        ImVec2(menu_width + 10.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(800.0f, 50.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(200.0f, 30.0f), ImVec2(800.0f, -1.0f));
    bool _message_menu_visible = true;

    ImGui::Begin("Message", &_message_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    static ImVec4 color_ok(0.0, 1.0, 0.0, 1.0);
    static ImVec4 color_error(1.0, 0.0, 0.0, 1.0);

    ImGui::TextColored(last_action_success ? color_ok : color_error, "%s",
        last_action_message.c_str());
    ImGui::End();

    // ---------------------------------------------------------------------------------
    float legends_width = 250.f * menu_scaling();
    ImGui::SetNextWindowPos(
        ImVec2(ImGui::GetIO().DisplaySize.x - legends_width, 0.0f),
        ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(200.0f, 30.0f), ImVec2(legends_width, -1.0f));
    bool _colors_menu_visible = true;

    ImGui::Begin("UI", &_colors_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    draw_legends();
    ImGui::End();

    // ------------------------------------------------------------------------
    ImGui::SetNextWindowPos(
        ImVec2(
            ImGui::GetIO().DisplaySize.x - legends_width - menu_width - 10, 0),
        ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _line_stack_menu_visible = true;

    ImGui::Begin("Line Stack", &_line_stack_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);
    draw_line_stack();
    ImGui::PopItemWidth();
    ImGui::End();
}

// //////////////////////////////////////////////////////////////////////////
// LEGENDS MENU
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_legends()
{
    float slider_width = ImGui::GetWindowWidth() * 0.25f;
    // ------------------------------------------------------------------
    // EDGES
    // ------------------------------------------------------------------
    {
        // Cast show_overlay to a bool because it is a unsigned int
        bool tmp = bool(viewer->data_list[edges_data_id].show_overlay);
        ImGui::Checkbox("##edges-UI", &tmp);
        viewer->data_list[edges_data_id].show_overlay = unsigned(tmp);

        ImGui::SameLine();
        if (double_color_edit("edge##UI", color_edge)) {
            recolor_edges();
        }
        ImGui::SameLine();
        ImGui::PushItemWidth(slider_width);
        float point_size
            = viewer->data_list[edges_data_id].point_size / pixel_ratio();
        if (ImGui::SliderFloat(
                "##edges_scaling", &point_size, 0.00f, 10.0f, "%1.f")) {
            viewer->data_list[edges_data_id].point_size
                = point_size * pixel_ratio();
        }
        ImGui::PopItemWidth();
        ImGui::SameLine();
        ImGui::Checkbox("vertex-id##edges-UI",
            &viewer->data_list[edges_data_id].show_vertid);
    }
    // ------------------------------------------------------------------
    // DISPL
    // ------------------------------------------------------------------
    {
        // Cast show_overlay to a bool because it is a unsigned int
        bool tmp = bool(viewer->data_list[displ_data_id].show_overlay);
        ImGui::Checkbox("##displ-UI", &tmp);
        viewer->data_list[displ_data_id].show_overlay = unsigned(tmp);

        ImGui::SameLine();
        if (double_color_edit("displ##UI", color_displ)) {
            recolor_displacements();
        }
        ImGui::SameLine();
        ImGui::PushItemWidth(slider_width);
        float point_size
            = viewer->data_list[displ_data_id].point_size / pixel_ratio();
        if (ImGui::SliderFloat(
                "##displ_point_scaling", &point_size, 0.00f, 10.0f, "%1.f")) {
            viewer->data_list[displ_data_id].point_size
                = point_size * pixel_ratio();
        }
        ImGui::PopItemWidth();
    }
    // ------------------------------------------------------------------
    // GRAD
    // ------------------------------------------------------------------
    {
        // Cast show_overlay to a bool because it is a unsigned int
        bool tmp = bool(viewer->data_list[gradient_data_id].show_overlay);
        ImGui::Checkbox("##grad-UI", &tmp);
        viewer->data_list[gradient_data_id].show_overlay = unsigned(tmp);

        ImGui::SameLine();
        if (double_color_edit("grad##UI", color_grad)) {
            recolor_grad_volume();
        }
        ImGui::SameLine();
        ImGui::PushItemWidth(slider_width);
        if (ImGui::SliderFloat(
                "##grad_scaling", &state.grad_scaling, 0.00f, 5.0f, "%1.f")) {
            redraw_grad_volume(state.use_opt_gradient);
        }
        ImGui::PopItemWidth();
    }

    // ------------------------------------------------------------------
    // OPT DISPL
    // ------------------------------------------------------------------
    {
        // Cast show_overlay to a bool because it is a unsigned int
        bool tmp = bool(viewer->data_list[opt_displ_data_id].show_overlay);
        ImGui::Checkbox("##opt displ-UI", &tmp);
        viewer->data_list[opt_displ_data_id].show_overlay = unsigned(tmp);

        ImGui::SameLine();
        if (double_color_edit("opt displ##UI", color_opt_displ)) {
            recolor_opt_displacements();
        }
    }
    ImGui::Separator();
    if (double_color_edit("selection##UI", color_sl)) {
        recolor_edges();
        recolor_displacements();
    };
}

// //////////////////////////////////////////////////////////////////////////
// IO MENU
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_io()
{
    // Scene
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen)) {
        float w = ImGui::GetContentRegionAvailWidth();
        float p = ImGui::GetStyle().FramePadding.x;
        if (ImGui::Button("Load##Scene", ImVec2((w - p) / 2.f, 0))) {
            load_scene();
        }
        ImGui::SameLine(0, p);
        if (ImGui::Button("Save##Scene", ImVec2((w - p) / 2.f, 0))) {
            save_scene();
        }
        ImGui::Checkbox(
            "convert to rigid bodies", &state.convert_to_rigid_bodies);
    }
}

// //////////////////////////////////////////////////////////////////////////
// EDIT MODE
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_edit_modes()
{
    if (ImGui::CollapsingHeader("Edit Mode", ImGuiTreeNodeFlags_DefaultOpen)) {
        float w = ImGui::GetContentRegionAvailWidth();
        float p = ImGui::GetStyle().FramePadding.x;

        for (uint i = 0; i < ViewerEditModeAll.size(); ++i) {
            bool needs_pop = false;
            if (edit_mode == ViewerEditModeAll[i]) {
                color_button();
                needs_pop = true;
            }
            char buf[100];
            sprintf(buf, "%s##Edit", ViewerEditModeNames[i].c_str());
            if (ImGui::Button(buf, ImVec2((w - p) / 2.f, 0))) {
                if (edit_mode == ViewerEditModeAll[i]) {
                    edit_mode = none;
                } else {
                    edit_mode = ViewerEditModeAll[i];
                }
            }
            if (needs_pop)
                ImGui::PopStyleColor(3);
            if (i % 2 == 0) {
                ImGui::SameLine(0, p);
            }
        }

        if (ImGui::Button("Connect##Edit", ImVec2(-1, 0))) {
            connect_selected_vertices();
        }

        if (ImGui::Button("Select All Vertices##Edit", ImVec2(-1, 0))) {
            state.selected_displacements.clear();
            Eigen::VectorXi all_idxs = Eigen::VectorXi::LinSpaced(
                state.vertices.rows(), 0, state.vertices.rows() - 1);
            state.selected_points = std::vector<int>(
                all_idxs.data(), all_idxs.data() + all_idxs.size());
            recolor_edges();
            recolor_displacements();
        }
        if (ImGui::Button("Select All Displacements##Edit", ImVec2(-1, 0))) {
            state.selected_points.clear();
            Eigen::VectorXi all_idxs = Eigen::VectorXi::LinSpaced(
                state.displacements.rows(), 0, state.displacements.rows() - 1);
            state.selected_displacements = std::vector<int>(
                all_idxs.data(), all_idxs.data() + all_idxs.size());
            recolor_edges();
            recolor_displacements();
        }

        if (ImGui::Button("Subdivide Edges##Edit", ImVec2(-1, 0))) {
            // Loop over edges creating a new vertex in the center
            const long num_old_vertices = state.vertices.rows();
            const long num_new_vertices = state.edges.rows();
            state.vertices.conservativeResize(
                num_old_vertices + num_new_vertices, state.vertices.cols());
            state.displacements.conservativeResize(
                num_old_vertices + num_new_vertices,
                state.displacements.cols());
            state.edges.conservativeResize(2 * state.edges.rows(), 2);

            auto fixed_dof = state.opt_problem.fixed_dof;
            fixed_dof.resize(num_old_vertices, 2);
            fixed_dof.conservativeResize(state.displacements.rows(), 2);

            for (long i = 0; i < num_new_vertices; i++) {
                state.vertices.row(num_old_vertices + i) = 0.5
                    * (state.vertices.row(state.edges(i, 1))
                        + state.vertices.row(state.edges(i, 0)));

                state.displacements.row(num_old_vertices + i) = 0.5
                    * (state.displacements.row(state.edges(i, 1))
                        + state.displacements.row(state.edges(i, 0)));

                for (int j = 0; j < 2; j++) {
                    fixed_dof(num_old_vertices + i, j)
                        = fixed_dof(state.edges(i, 0), j)
                        && fixed_dof(state.edges(i, 1), j);
                }

                state.edges.row(num_new_vertices + i) = state.edges.row(i);
                state.edges(i, 1) = num_old_vertices + i;
                state.edges(num_new_vertices + i, 0) = num_old_vertices + i;
            }

            fixed_dof.resize(fixed_dof.size(), 1);

            state.reset_scene();
            state.opt_problem.fixed_dof = fixed_dof;

            state_history.push_back(state);
            load_state();
        }

        if (ImGui::Button("Smooth Vertices##Edit", ImVec2(-1, 0))) {
            // Loop over edges creating a new vertex in the center
            double w0 = 0.75; // TODO: Expose this parameter
            Eigen::MatrixX2d smoothed_vertices = w0 * state.vertices;
            Eigen::MatrixX2d smoothed_displacements = w0 * state.displacements;

            auto adjacency_list = state.create_adjacency_list();

            for (long i = 0; i < state.vertices.rows(); i++) {
                const unsigned long num_neighbors = adjacency_list[i].size();
                for (const auto& vertex_idx : adjacency_list[i]) {
                    smoothed_vertices.row(i) += ((1 - w0) / num_neighbors)
                        * state.vertices.row(vertex_idx);
                    smoothed_displacements.row(i) += ((1 - w0) / num_neighbors)
                        * state.displacements.row(vertex_idx);
                }
            }

            state.vertices = smoothed_vertices;
            state.displacements = smoothed_displacements;

            auto fixed_dof = state.opt_problem.fixed_dof;
            state.reset_scene();
            state.opt_problem.fixed_dof = fixed_dof;

            state_history.push_back(state);
            load_state();
        }

        if (ImGui::Button("Remove Free Vertices##Edit", ImVec2(-1, 0))) {
            if (state.remove_free_vertices()) {
                state_history.push_back(state);
                load_state();
            }
        }

        if ((state.selected_points.size() > 0
                || state.selected_displacements.size() > 0)
            && ImGui::Button("Select Connected##Edit", ImVec2(-1, 0))) {
            std::vector<int>& selection = state.selected_points.size() > 0
                ? state.selected_points
                : state.selected_displacements;
            // Build adjacency list
            std::vector<std::vector<int>> adjacency_list
                = state.create_adjacency_list();

            std::unordered_set<int> new_selection;
            for (const auto& selected_id : selection) {
                state.find_connected_vertices(
                    selected_id, adjacency_list, new_selection);
            }
            selection.assign(new_selection.begin(), new_selection.end());

            recolor_edges();
            recolor_displacements();
        }

        if (state.selected_points.size() > 0
            && ImGui::Button("Duplicate Selected##Edit", ImVec2(-1, 0))) {
            // Duplicate selected vertices and edges that have both end-points
            // selected
            long num_old_vertices = state.vertices.rows();
            Eigen::Vector2d delta_com;
            delta_com << 0, -1;
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
    }

    // Menu for fixing vertex positions
    if (state.selected_points.size() > 0
        && ImGui::CollapsingHeader(
            "Static Vertices##static", ImGuiTreeNodeFlags_DefaultOpen)) {
        // Initial button state is all(fixed_dof(selected_points))
        bool x_fixed_originally = true, y_fixed_originally = true;
        for (int point : state.selected_points) {
            x_fixed_originally &= state.opt_problem.fixed_dof(point);
            y_fixed_originally &= state.opt_problem.fixed_dof(
                point + state.displacements.rows());
        }
        bool x_fixed = x_fixed_originally, y_fixed = y_fixed_originally;

        ImGui::Checkbox("fixed x position##static", &x_fixed);
        if (x_fixed != x_fixed_originally) {
            for (int point : state.selected_points) {
                state.opt_problem.fixed_dof(point) = x_fixed;
                if (x_fixed) {
                    state.displacements(point, 0) = 0.0;
                }
            }
            state_history.push_back(state);
            load_state();
            recolor_edges();
        }

        ImGui::Checkbox("fixed y position##static", &y_fixed);
        if (y_fixed != y_fixed_originally) {
            for (int point : state.selected_points) {
                state.opt_problem.fixed_dof(point + state.displacements.rows())
                    = y_fixed;
                if (y_fixed) {
                    state.displacements(point, 1) = 0.0;
                }
            }
            state_history.push_back(state);
            load_state();
            recolor_edges();
        }
    }
}

void ViewerMenu::draw_line_stack()
{
    static int num_lines = 3;
    static double scale_displacment = 10;
    ImGui::InputIntBounded("line count##line-stack", &num_lines, 0,
        std::numeric_limits<int>::max(), 1, 10);
    ImGui::InputDouble("scale disp.##line-stack", &scale_displacment);
    if (ImGui::Button("Make Line Stack##Edit", ImVec2(-1, 0))) {
        Eigen::MatrixX2d vertices(2 * num_lines + 2, 2);
        Eigen::MatrixX2d displacements
            = Eigen::MatrixX2d::Zero(2 * num_lines + 2, 2);
        Eigen::MatrixX2i edges(num_lines + 1, 2);

        vertices.row(0) << -0.05, 0.1;
        vertices.row(1) << 0.05, 0.2;
        displacements(0, 1) = -1;
        displacements(1, 1) = -1;
        displacements *= scale_displacment;
        edges.row(0) << 0, 1;

        Eigen::VectorXd ys = Eigen::VectorXd::LinSpaced(num_lines, 0, -1);
        for (int i = 1; i < edges.rows(); i++) {
            vertices.row(2 * i) << -1, ys(i - 1);
            vertices.row(2 * i + 1) << 1, ys(i - 1);
            edges.row(i) << 2 * i, 2 * i + 1;
        }

        state.vertices = vertices;
        state.displacements = displacements;
        state.edges = edges;
        state.reset_scene();
        state_history.push_back(state);
        load_state();
    }
}

// //////////////////////////////////////////////////////////////////////////
// CCD STEPS
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_ccd_steps()
{
    if (ImGui::CollapsingHeader("CCD", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::SliderFloat("time", &(state.current_time), 0.0, 1.0)) {
            redraw_at_time();
        }

        if (ImGui::Button("Run CCD", ImVec2(-1, 0))) {
            compute_collisions();
        }
        if (state.volumes.size() > 0) {
            if (ImGui::InputInt("volume##volume", &state.current_volume)) {
                redraw_grad_volume(/*opt_gradient=*/false);
            }
        }
        //        ImGui::Checkbox("skip empty", &state.skip_no_impact_edge);
        ImGui::Text("||jac_j(i)|| = \t%.3g", state.get_volume_grad().norm());
    }
}

void ViewerMenu::draw_collision_menu()
{

    int idx_ctr_type = static_cast<int>(state.constraint_function);

    if (ImGui::Combo("Constraint##opt", &idx_ctr_type, ccd::ConstraintNames,
            CCD_IM_ARRAYSIZE(ccd::ConstraintNames))) {
        state.constraint_function
            = static_cast<ccd::ConstraintType>(idx_ctr_type);
    }
    {
        ImGui::Indent();
        ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);
        if (state.constraint_function == ccd::ConstraintType::VOLUME) {
            volume_constraint_menu(state.volume_constraint);
        } else {
            barrier_constraint_menu(state.barrier_constraint);
        }
        ImGui::PopItemWidth();
        ImGui::Unindent();
    }
}

// //////////////////////////////////////////////////////////////////////////
// OPTIMIZATION
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_optimization()
{
    using namespace opt;
    int idx_optimization_method = static_cast<int>(state.opt_method);

    if (ImGui::CollapsingHeader(
            "Collision Optimization", ImGuiTreeNodeFlags_DefaultOpen)) {

        draw_collision_menu();

        if (ImGui::Combo("method##opt", &idx_optimization_method,
                ccd::OptimizationMethodNames,
                CCD_IM_ARRAYSIZE(ccd::OptimizationMethodNames))) {
            state.opt_method
                = static_cast<ccd::OptimizationMethod>(idx_optimization_method);
        }

        {
            ImGui::Indent();
            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);
            switch (state.opt_method) {
            case ccd::OptimizationMethod::LINEARIZED_CONSTRAINTS:
                solver_menu(state.qp_solver);
                break;
            case ccd::OptimizationMethod::NCP:
                solver_menu(state.ncp_solver);
                break;
            case ccd::OptimizationMethod::IPOPT:
                solver_menu(state.ipopt_solver);
                break;
            case ccd::OptimizationMethod::NLOPT:
                solver_menu(state.nlopt_solver);
                break;
            case ccd::OptimizationMethod::BARRIER_NEWTON:
                solver_menu(state.barrier_newton_solver);
                break;
            }
            ImGui::PopItemWidth();
            ImGui::Unindent();
        }

        ImGui::Checkbox(
            "continue optimization##opt", &(state.reuse_opt_displacements));

        if (ImGui::Button("Optimize##opt", ImVec2(-1, 0))) {
            optimize_displacements();
        }

        if (ImGui::Button("Process##opt", ImVec2(-1, 0))) {
            post_process_optimization();
        }
    }
} // namespace ccd

// //////////////////////////////////////////////////////////////////////////
// OPTIMIZATION RESULTS
// //////////////////////////////////////////////////////////////////////////
void ViewerMenu::draw_optimization_results()
{
    if (!ImGui::CollapsingHeader(
            "Opt Results", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }
    if (state.opt_results.finished) {
        // -------------------------------------------------------------------
        // DETAILS
        // -------------------------------------------------------------------

        ImGui::BeginChild("Opt Results Detail",
            ImVec2(ImGui::GetWindowContentRegionWidth() * 0.9f, 100), false);

        ImGui::Text("method = %s",
            ccd::OptimizationMethodNames[static_cast<int>(state.opt_method)]);

        double energy = state.get_opt_functional();
        if (energy > -1) {
            ImGui::Text("energy = %.3g", energy);
        } else {
            ImGui::Text("energy = NO_INFO");
        }

        ImGui::Text("displacements");
        Eigen::MatrixX2d disp = state.get_opt_displacements();
        ImGui::Columns(2, /*id=*/nullptr, /*border=*/false);
        for (uint i = 0; i < disp.rows(); i++) {
            ImGui::Text("%.3g", disp(i, 0));
            ImGui::NextColumn();
            ImGui::Text("%.3g", disp(i, 1));
            ImGui::NextColumn();
        }
        ImGui::Columns(1);
        ImGui::EndChild();

        // -------------------------------------------------------------------
        // EXPLORE
        // -------------------------------------------------------------------
        ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.7f);
        if (ImGui::SliderFloat(
                "time##opt-results", &(state.current_opt_time), 0.0, 1.0)) {
            redraw_at_opt_time();
        }
        ImGui::PopItemWidth();
        if (state.u_history.size() > 0) {
            if (ImGui::InputIntBounded("step##opt-results",
                    &(state.current_opt_iteration), -1,
                    int(state.u_history.size() - 1), 1, 10)) {

                redraw_opt_displacements();
                recolor_opt_collisions();
                redraw_grad_volume(/*opt_gradient=*/true);
                redraw_at_opt_time();
            }
        }
        if (state.g_history.size() > 0)
            if (ImGui::InputInt("volume##ee_opt", &state.current_volume)) {
                redraw_grad_volume(/*opt_gradient=*/true);
            }
    }

    // -------------------------------------------------------------------
    // LOAD/SAVE
    // -------------------------------------------------------------------
    float w = ImGui::GetContentRegionAvailWidth();
    float p = ImGui::GetStyle().FramePadding.x;
    if (ImGui::Button("Load##opt-results", ImVec2((w - p) / 2.f, 0))) {
        load_optimization();
    }
    ImGui::SameLine(0, p);
    if (ImGui::Button("Save##opt-results", ImVec2((w - p) / 2.f, 0))) {
        save_optimization();
    }
}

} // namespace ccd
