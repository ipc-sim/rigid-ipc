#include "UISimState.hpp"

#include <viewer/imgui_ext.hpp>
#include <logger.hpp>

#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <imgui.h>

namespace ipc::rigid {

void UISimState::draw_menu()
{
    static bool player_menu = true;
    static bool settings_menu = true;
    static bool collisions_menu = true;
    static bool legends_menu = true;

    //--------------------------------------------------------------------------
    {
        draw_io();
        if (m_has_scene) {
            if (ImGui::CollapsingHeader(
                    "Player", &player_menu, ImGuiTreeNodeFlags_DefaultOpen)) {
                draw_simulation_player();
            }

            ImGui::SetNextItemOpen(false, ImGuiCond_FirstUseEver);
            if (ImGui::CollapsingHeader("Settings", &settings_menu)) {
                draw_settings();
            }
        }
    }
} // namespace ipc::rigid

void UISimState::draw_io()
{
    ImGui::PushItemWidth(200);
    const char* level_strings[] = { "trace", "debug",    "info", "warning",
                                    "error", "critical", "off" };
    int log_level = spdlog::get_level();
    if (ImGui::Combo(
            "Log level##logger", &log_level, level_strings,
            CCD_IM_ARRAYSIZE(level_strings))) {
        set_logger_level(static_cast<spdlog::level::level_enum>(log_level));
    }
    ImGui::PopItemWidth();

    if (ImGui::Button("Load##IO")) {
        std::string fname = igl::file_dialog_open();
        if (fname != "") {
            load(fname);
        }
    }

    if (m_has_scene) {
        ImGui::SameLine();

        if (ImGui::Button("Reload##IO")) {
            reload();
        }

        ImGui::SameLine();
        if (ImGui::Button("Save checkpoint")) {
            std::string fname = igl::file_dialog_save();
            save(fname);
        }

        ImGui::SameLine();
        if (ImGui::Button("Save GLTF")) {
            std::string filename = igl::file_dialog_save();
            if (filename != "") {
                save_gltf(filename);
            }
        }

        ImGui::SameLine();

        if (m_is_gif_recording) {
            if (ImGui::Button("Save GIF")) {
                end_recording();
            }
        } else if (ImGui::Button("Record GIF")) {
            std::string fname = igl::file_dialog_save();
            start_recording(fname);
        }

        ImGui::Text("Scene file:");
        ImGui::BeginChild(
            "##filename", ImVec2(-1, 2.5 * ImGui::GetFontSize()), false,
            ImGuiWindowFlags_HorizontalScrollbar);
        {
            ImGui::Text("%s", m_state.scene_file.c_str());
        }
        ImGui::EndChild();
    }
}

void UISimState::draw_simulation_player()
{
    PlayerState player_state = m_player_state;

    if (ImGui::Button(
            player_state == PlayerState::Playing ? "Pause" : "Play")) {
        player_state = (m_player_state == PlayerState::Playing)
            ? PlayerState::Paused
            : PlayerState::Playing;
        if (player_state == PlayerState::Paused) {
            log_simulation_time();
        }
        if (player_state == PlayerState::Playing
            && m_state.m_num_simulation_steps
                >= m_state.m_max_simulation_steps) {
            m_state.m_max_simulation_steps = -1; // Turn off breaking
        }
    }
    m_player_state = static_cast<PlayerState>(player_state);

    ImGui::SameLine();
    if (ImGui::Button("Step##SimPlayer")) {
        m_player_state = PlayerState::Playing;
        pre_draw_loop();
        m_player_state = PlayerState::Paused;
    }

    ImGui::SameLine();

    ImGui::PushItemWidth(200);
    if (ImGui::SliderInt(
            "Frame##Replay", &m_state.m_num_simulation_steps, 0,
            m_state.state_sequence.size() - 1)) {
        m_state.problem_ptr->state(
            m_state.state_sequence[m_state.m_num_simulation_steps]);
        redraw_scene();
        m_scene_changed = true;
    }
    ImGui::PopItemWidth();

    // --------------------------------------------------------------------

    if (ImGui::TreeNode("Breakpoints")) {
        ImGui::Checkbox("Intersecting", &m_bkp_has_intersections);
        ImGui::SameLine();
        ImGui::HelpMarker("yes - stop playing if step has intersections.");

        ImGui::Checkbox("Optimization failed", &m_bkp_optimization_failed);
        ImGui::SameLine();
        ImGui::HelpMarker("yes - stop playing if step optimization failed.");

        ImGui::Checkbox("On collision", &m_bkp_had_collision);
        ImGui::SameLine();
        ImGui::HelpMarker("yes - stop playing if step had a collision.");

        ImGui::TreePop();
    }
}

void UISimState::draw_settings()
{
    auto config = m_state.get_active_config();
    ImGui::BeginChild(
        "##config", ImVec2(-1, 300), false,
        ImGuiWindowFlags_HorizontalScrollbar);
    {
        ImGui::TreeNodeJson(config);
    }
    ImGui::EndChild();
}

} // namespace ipc::rigid
