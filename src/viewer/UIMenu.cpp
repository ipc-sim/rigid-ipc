#include "UISimState.hpp"

#include <viewer/imgui_ext.hpp>
#include <viewer/solver_view.hpp>
#include <viewer/constraint_view.hpp>

#include <logger.hpp>

namespace ccd {

void UISimState::draw_menu()
{
    draw_labels_window();

    float menu_width = 200.f * menu_scaling();

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
            draw_simulation_player();
            draw_ccd();
        }
    }
    ImGui::PopItemWidth();
    ImGui::End();
}

void UISimState::draw_io()
{
    if (ImGui::Button("Load##IO", ImVec2(-1, 0))) {
        std::string fname = igl::file_dialog_open();
        load(fname);
    }
}

void UISimState::draw_simulation_player()
{
    int player_state = static_cast<int>(m_player_state);
    ImGui::InputDouble(
        "dt##SimPlayer", &m_state.m_timestep_size, 0.01, 0.1, "%.3g");
    ImGui::RadioButton("play##SimPlayer", &player_state, PlayerState::Playing);
    ImGui::RadioButton("pause##SimPlayer", &player_state, PlayerState::Paused);
    m_player_state = static_cast<PlayerState>(player_state);

    if (ImGui::Button("Step##SimPlayer", ImVec2(-1, 0))) {
        simulation_step();
    }

    ImGui::Checkbox("break had collision", &m_bkp_had_collision);
    ImGui::Checkbox("break has collision", &m_bkp_has_collision);
    ImGui::Text("Step %i", m_state.m_num_simulation_steps);
}

void UISimState::draw_ccd() {
    if(ImGui::CollapsingHeader("Solver")){
        solver_menu(m_state.m_ccd_solver);
    }
    if(ImGui::CollapsingHeader("Constraint")){
        barrier_constraint_menu(m_state.m_ccd_constraint);
    }
}
} // namespace ccd
