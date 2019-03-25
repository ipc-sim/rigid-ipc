#include "viewer.hpp"

#define CCD_IM_ARRAYSIZE(_ARR) (int(sizeof(_ARR) / sizeof(*_ARR)))

namespace ccd {

void ViewerMenu::draw_menu()
{
    float menu_width = 180.f * menu_scaling();
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
    draw_edit_modes();
    draw_ui_settings();
    draw_ccd_steps();

    ImGui::PopItemWidth();
    ImGui::End();

    // ---------------------------------------------------------------------------------
    ImGui::SetNextWindowPos(
        ImVec2(menu_width + 10.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _opt_menu_visible = true;

    ImGui::Begin("Optimization", &_opt_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);

    draw_optimization();
    draw_optimization_results();

    ImGui::PopItemWidth();
    ImGui::End();

    // ---------------------------------------------------------------------------------
    ImGui::SetNextWindowPos(
        ImVec2(2 * menu_width + 10.0f, 0.0f), ImGuiSetCond_FirstUseEver);
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
}

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
    }
}

void color_button()
{
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5, 0.0, 0.0, 1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f, 0.0, 0.0, 1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(1.0, 0.0, 0.0, 1.0));
}

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
                edit_mode = ViewerEditModeAll[i];
            }
            if (needs_pop)
                ImGui::PopStyleColor(3);
            if (i % 2 == 0) {
                ImGui::SameLine(0, p);
            }
        }
    }
    if (ImGui::Button("Connect##Edit", ImVec2(-1, 0))) {
        connect_selected_vertices();
    }
}

void ViewerMenu::draw_ui_settings()
{
    if (ImGui::CollapsingHeader(
            "UI Settings", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Checkbox(
            "show surface", &(viewer->data_list[surface_data_id].show_overlay));

        ImGui::Checkbox("show vertex id",
            &(viewer->data_list[surface_data_id].show_vertid));
        ImGui::InputFloat("point size ##surface",
            &(viewer->data_list[surface_data_id].point_size));

        ImGui::Checkbox("show displacements",
            &(viewer->data_list[displ_data_id].show_overlay));
        ImGui::InputFloat("point size ##displ",
            &(viewer->data_list[displ_data_id].point_size));

        ImGui::Checkbox(
            "show volumes", &(viewer->data_list[volume_data_id].show_faces));
        ImGui::Checkbox("show volumes border",
            &(viewer->data_list[volume_data_id].show_lines));

        ImGui::Checkbox("show gradient",
            &(viewer->data_list[gradient_data_id].show_overlay));

        if (ImGui::SliderFloat("time", &(state.time), 0.0, 1.0)) {
            redraw_at_time();
        }
    }
}

void ViewerMenu::draw_ccd_steps()
{
    static int idx_detection_method = 0;
    if (ImGui::CollapsingHeader("CCD Steps", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Combo("method##ccd", &idx_detection_method,
                ccd::DetectionMethodNames,
                CCD_IM_ARRAYSIZE(ccd::DetectionMethodNames))) {
            state.detection_method
                = static_cast<ccd::DetectionMethod>(idx_detection_method);
        }
        ImGui::InputDouble("vol. epsilon", &state.volume_epsilon);

        if (ImGui::Button("Detect EV Collisions", ImVec2(-1, 0))) {
            detect_edge_vertex_collisions();
        }

        if (state.ev_impacts.size()) {
            if (ImGui::InputInt(
                    "VE Impact##ev_impact", &state.current_ev_impact)) {
                goto_ev_impact(state.current_ev_impact);
            }
        }
        if (state.ee_impacts.size()) {
            if (ImGui::InputInt("EE Impact##ee_impact", &state.current_edge)) {
                goto_ee_impact(state.current_edge);
            }
            ImGui::Checkbox("skip empty", &state.skip_no_impact_edge);
        }
        ImGui::Separator();

        if (state.ev_impacts.size()) {
            if (ImGui::Button("Compute Collision Volumes", ImVec2(-1, 0))) {
                compute_collision_volumes();
            }
        }
    }
}

void ViewerMenu::draw_optimization()
{
    int idx_optimization_method = state.solver_settings.method;
    int idx_qp_solver = state.solver_settings.qp_solver;

    if (ImGui::CollapsingHeader(
            "Displacement Optimization", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Combo("method##opt", &idx_optimization_method,
                ccd::opt::OptimizationMethodNames,
                CCD_IM_ARRAYSIZE(ccd::opt::OptimizationMethodNames))) {
            state.solver_settings.method
                = static_cast<ccd::opt::OptimizationMethod>(
                    idx_optimization_method);
        }

        if (state.solver_settings.method == opt::LINEARIZED_CONSTRAINTS) {

            if (ImGui::Combo("QP solver##opt", &idx_qp_solver,
                    ccd::opt::QPSolverNames,
                    CCD_IM_ARRAYSIZE(ccd::opt::QPSolverNames))) {
                state.solver_settings.qp_solver
                    = static_cast<ccd::opt::QPSolver>(idx_qp_solver);
            }
        }

        ImGui::InputInt("max iter##opt", &state.solver_settings.max_iter);
        ImGui::InputInt("verbosity##opt", &state.solver_settings.verbosity);

        ImGui::Checkbox(
            "continue optimization##opt", &(state.reuse_opt_displacements));

        ImGui::Checkbox(
            "recompute col. set##opt", &(state.recompute_collision_set));

        if (ImGui::Button("Optimize##opt", ImVec2(-1, 0))) {
            optimize_displacements();
            update_graph(
                surface_data_id, state.get_opt_vertex_at_time(-1), state.edges);
        }
    }
}

void ViewerMenu::draw_optimization_results()
{
    if (ImGui::CollapsingHeader(
            "Opt Results", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (state.opt_results.finished) {
            ImGui::BeginChild("Opt Results Detail",
                ImVec2(ImGui::GetWindowContentRegionWidth() * 0.9f, 100),
                false);
            ImGui::Text("method = %s",
                ccd::opt::OptimizationMethodNames[state.opt_results.method]);
            ImGui::Text(
                "energy = %.3g", state.get_opt_functional(state.opt_iteration));

            ImGui::Text("displacements");
            Eigen::MatrixX2d disp
                = state.get_opt_displacements(state.opt_iteration);
            ImGui::Columns(2, /*id=*/nullptr, /*border=*/false);
            for (uint i = 0; i < disp.rows(); i++) {
                ImGui::Text("%.3g", disp(i, 0));
                ImGui::NextColumn();
                ImGui::Text("%.3g", disp(i, 1));
                ImGui::NextColumn();
            }
            ImGui::Columns(1);
            ImGui::EndChild();

            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.7f);
            if (ImGui::SliderFloat(
                    "time##opt-results", &(state.opt_time), 0.0, 1.0)) {
                redraw_at_opt_time();
            }
            ImGui::PopItemWidth();
            if (state.u_history.size() > 0) {
                if (ImGui::InputInt(
                        "step##opt-results", &(state.opt_iteration), 1, 10)) {
                    redraw_opt_displacements();
                    redraw_at_opt_time();
                }
            }
        }

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
}

} // namespace ccd
