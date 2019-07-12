#include "UISimState.hpp"

#include <viewer/imgui_ext.hpp>

#include <logger.hpp>

namespace ccd {

void UISimState::draw_menu()
{
    draw_labels_window();

    float menu_width = 200.f * menu_scaling();
    static bool player_menu = true;
    static bool settings_menu = true;
    static bool collisions_menu = true;

    //--------------------------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _viewer_menu_visible = true;

    ImGui::Begin("Viewer", &_viewer_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);
    {

        const char* level_strings[] = SPDLOG_LEVEL_NAMES;
        if (ImGui::Combo("log-level##logger", &m_log_level, level_strings,
                CCD_IM_ARRAYSIZE(level_strings))) {
            spdlog::set_level(
                static_cast<spdlog::level::level_enum>(m_log_level));
        }

        draw_io();
        if (m_has_scene) {

            if (ImGui::CollapsingHeader(
                    "Player", &player_menu, ImGuiTreeNodeFlags_DefaultOpen)) {
                draw_simulation_player();
            }
            if (ImGui::CollapsingHeader("Collisions", &collisions_menu,
                    ImGuiTreeNodeFlags_DefaultOpen)) {
                draw_collision_menu();
            }
            if (ImGui::CollapsingHeader(
                    "Settings", &settings_menu, ImGuiTreeNodeFlags_DefaultOpen))
                draw_settings();
        }
    }
    ImGui::PopItemWidth();
    ImGui::End();

    //--------------------------------------------------------------------------
    if (m_has_scene) {
        float legends_width = 250.f * menu_scaling();
        ImGui::SetNextWindowPos(
            ImVec2(ImGui::GetIO().DisplaySize.x - legends_width, 0.0f),
            ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints(
            ImVec2(200.0f, 30.0f), ImVec2(legends_width, -1.0f));
        bool _colors_menu_visible = true;

        ImGui::Begin("UI", &_colors_menu_visible,
            ImGuiWindowFlags_NoSavedSettings
                | ImGuiWindowFlags_AlwaysAutoResize);
        {
            draw_legends();
            if (ImGui::Checkbox("show as position-deltas", &m_show_as_delta)) {
                redraw_scene();
            }
        }
        ImGui::End();
    }
} // namespace ccd

void UISimState::draw_io()
{
    if (ImGui::Button("Load##IO", ImVec2(-1, 0))) {
        std::string fname = igl::file_dialog_open();
        load(fname);
    }
    if (m_has_scene) {
        if (ImGui::Button("Reload##IO", ImVec2(-1, 0))) {
            reload();
        }
        ImGui::BeginChild("##filename",
            ImVec2(ImGui::GetWindowContentRegionWidth() * 0.9f,
                ImGui::GetFontSize() * 3),
            false, ImGuiWindowFlags_HorizontalScrollbar);
        {
            ImGui::Text("%s", m_state.scene_file.c_str());
            ImGui::SetScrollHere(1.0f);
        }
        ImGui::EndChild();
    }
}

void UISimState::draw_simulation_player()
{
    int player_state = static_cast<int>(m_player_state);

    ImGui::RadioButton("play##SimPlayer", &player_state, PlayerState::Playing);
    ImGui::RadioButton("pause##SimPlayer", &player_state, PlayerState::Paused);
    m_player_state = static_cast<PlayerState>(player_state);

    if (ImGui::Button("Step##SimPlayer", ImVec2(-1, 0))) {
        simulation_step();
    }
    // --------------------------------------------------------------------
    ImGui::Checkbox("auto. solve collisions", &m_state.m_solve_collisions);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - solve collisions automatically on each step.");

    ImGui::Checkbox("break had collision", &m_bkp_had_collision);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - stop playing if step had a collision.");

    ImGui::Checkbox("break has collision", &m_bkp_has_collision);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - stop playing if step has unsolved collision.");

    // --------------------------------------------------------------------
    ImGui::Text("Step %i", m_state.m_num_simulation_steps);

    int item_current = m_show_next_step ? 1 : 0;
    if (ImGui::Combo(
            "##prev-next", &item_current, "previous_step\0next_step\0\0")) {
        m_show_next_step = item_current == 1 ? true : false;
        redraw_scene();
    }
    if (ImGui::DragDouble("interval time", &m_interval_time, 0.1, 0.0, 1.0)) {
        redraw_scene();
    }

    ImGui::SameLine();
    ImGui::HelpMarker("yes - shows velocities and forces as position-delta\n "
                      "i.e scaling them by time.");
}

void UISimState::draw_collision_menu()
{
    if (ImGui::Button("Solve Collisions", ImVec2(-1, 0))) {
        solve_collisions();
    }
    if (ImGui::Button("Step Solve Col. ", ImVec2(-1, 0))) {
        step_solve_collisions();
    }
    if (m_state.problem_ptr->has_barrier_constraint()) {
        double eps = m_state.problem_ptr->get_barrier_epsilon();
        if (ImGui::InputDouble("epsilon", &eps, 1e-3, 0.1, "%.3g")) {
            m_state.problem_ptr->set_barrier_epsilon(eps);
        }
    }
    ImGui::Text(
        "Outer it.: %i", m_state.ccd_solver_ptr->num_outer_iterations());
}
void UISimState::draw_settings()
{
    ImGui::Text("timestep_size: %.3g", m_state.m_timestep_size);
    Eigen::VectorXd G = m_state.problem_ptr->gravity();
    if (G.rows() == 2) {
        ImGui::Text("gravity: %.3g, %.3g", G(0), G(1));
    } else {
        ImGui::Text("gravity: %.3g, %.3g, %.3g", G(0), G(1), G(2));
    }
}

void UISimState::draw_legends()
{
    float slider_width = ImGui::GetWindowWidth() * 0.25f;

    for (auto& label : get_data_names()) {
        auto ptr = get_data(label);

        bool tmp = bool(ptr->data().show_overlay);
        ImGui::Checkbox(("##UI-" + label).c_str(), &tmp);
        ptr->data().show_overlay = tmp;
        ImGui::SameLine();
        if (ImGui::DoubleColorEdit3((label + "##UI").c_str(), ptr->m_color)) {
            ptr->recolor();
        }
        if (ptr->is_graph()) {
            ImGui::SameLine();
            ImGui::PushItemWidth(slider_width);
            {
                float point_size = ptr->data().point_size / pixel_ratio();
                if (ImGui::SliderFloat(("##UI-scaling" + label).c_str(),
                        &point_size, 0.00f, 10.0f, "%1.f")) {
                    ptr->data().point_size = point_size * pixel_ratio();
                }
            }
            ImGui::PopItemWidth();
            ImGui::SameLine();
            if (ImGui::Checkbox(
                    ("data##UI-" + label).c_str(), &ptr->show_vertex_data)) {
                ptr->update_vertex_data();
            }
        }
    }
}

} // namespace ccd
