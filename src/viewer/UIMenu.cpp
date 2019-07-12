#include "UISimState.hpp"

#include <viewer/imgui_ext.hpp>

#include <logger.hpp>

namespace ccd {

void UISimState::draw_menu()
{
    draw_labels_window();

    float menu_width = 200.f * menu_scaling();

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
            draw_simulation_player();
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
        }
        ImGui::End();
    }
}

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
    ImGui::Checkbox("break had collision", &m_bkp_had_collision);
    ImGui::Checkbox("break has collision", &m_bkp_has_collision);
    ImGui::Text("Step %i", m_state.m_num_simulation_steps);
    if (ImGui::DragDouble("t", &m_interval_time, 0.1, 0.0, 1.0)) {
        redraw_scene();
    }
    ImGui::Separator();
    if (ImGui::Checkbox("vel and Fc as deltas", &m_show_as_delta)){
        redraw_scene();
    }
    ImGui::Separator();
    ImGui::Text("timestep_size: %.3g", m_state.m_timestep_size);
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
            if (ImGui::Checkbox(("data##UI-"+ label).c_str(), &ptr->show_vertex_data)){
                ptr->update_vertex_data();
            }
        }
    }
}

} // namespace ccd
