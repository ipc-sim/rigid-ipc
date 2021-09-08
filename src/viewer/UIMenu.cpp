#include "UISimState.hpp"

#include <igl/Timer.h>

#include <viewer/imgui_ext.hpp>

#include <logger.hpp>

namespace ipc::rigid {

void UISimState::draw_menu()
{
    // if (m_show_vertex_data) {
    //     draw_labels_window();
    // }

    float menu_width = 220.f * menu_scaling();
    static bool player_menu = true;
    static bool settings_menu = true;
    static bool collisions_menu = true;

    //--------------------------------------------------------------------------
    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(
        ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _viewer_menu_visible = true;

    ImGui::Begin(
        "Controls", &_viewer_menu_visible,
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);
    {
        const char* level_strings[] = { "trace", "debug",    "info", "warning",
                                        "error", "critical", "off" };
        int log_level = spdlog::get_level();
        if (ImGui::Combo(
                "log-level##logger", &log_level, level_strings,
                CCD_IM_ARRAYSIZE(level_strings))) {
            set_logger_level(static_cast<spdlog::level::level_enum>(log_level));
        }

        draw_io();
        if (m_has_scene) {

            if (ImGui::CollapsingHeader(
                    "Player", &player_menu, ImGuiTreeNodeFlags_DefaultOpen)) {
                draw_simulation_player();
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
        float legends_width = 270.f * menu_scaling();
        ImGui::SetNextWindowPos(
            ImVec2(ImGui::GetIO().DisplaySize.x - legends_width, 0.0f),
            ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints(
            ImVec2(200.0f, 30.0f), ImVec2(legends_width, -1.0f));
        bool _colors_menu_visible = true;

        ImGui::Begin(
            "Legend", &_colors_menu_visible,
            ImGuiWindowFlags_NoSavedSettings
                | ImGuiWindowFlags_AlwaysAutoResize);
        {
            draw_legends();
        }
        ImGui::End();
    }
} // namespace ipc::rigid

void UISimState::draw_io()
{
    if (ImGui::Button("Load##IO", ImVec2(-1, 0))) {
        std::string fname = igl::file_dialog_open();
        if (fname != "") {
            load(fname);
        }
    }
    if (m_has_scene) {
        if (ImGui::Button("Reload##IO", ImVec2(-1, 0))) {
            reload();
        }
        ImGui::BeginChild(
            "##filename",
            ImVec2(
                ImGui::GetWindowContentRegionWidth() * 0.9f,
                ImGui::GetFontSize() * 3),
            false, ImGuiWindowFlags_HorizontalScrollbar);
        {
            ImGui::Text("%s", m_state.scene_file.c_str());
            ImGui::SetScrollHere(1.0f);
        }
        ImGui::EndChild();

        if (ImGui::Button("Save simulation", ImVec2(-1, 0))) {
            std::string fname = igl::file_dialog_save();
            save(fname);
        }

        // if (ImGui::Button("Save OBJ sequence", ImVec2(-1, 0))) {
        //     std::string dir_name = igl::file_dialog_save();
        //     if (dir_name != "") {
        //         save_obj_sequence(dir_name);
        //     }
        // }

        if (ImGui::Button("Save gltf", ImVec2(-1, 0))) {
            std::string filename = igl::file_dialog_save();
            if (filename != "") {
                save_gltf(filename);
            }
        }

        if (ImGui::Button("Save screenshot", ImVec2(-1, 0))) {
            std::string fname = igl::file_dialog_save();
            save_screenshot(fname);
        }
        if (m_is_gif_recording) {
            if (ImGui::Button("End recording", ImVec2(-1, 0))) {
                end_recording();
            }
        } else if (ImGui::Button("Start recording", ImVec2(-1, 0))) {
            std::string fname = igl::file_dialog_save();
            start_recording(fname);
        }
    }
}

void UISimState::draw_simulation_player()
{
    int player_state = static_cast<int>(m_player_state);

    ImGui::RadioButton("play##SimPlayer", &player_state, PlayerState::Playing);
    ImGui::RadioButton("pause##SimPlayer", &player_state, PlayerState::Paused);
    if (m_player_state == PlayerState::Playing
        && player_state == PlayerState::Paused && !replaying) {
        log_simulation_time();
    }
    if (m_state.m_num_simulation_steps >= m_state.m_max_simulation_steps
        && m_player_state == PlayerState::Paused
        && player_state == PlayerState::Playing) {
        m_state.m_max_simulation_steps = -1; // Turn off breaking
    }
    m_player_state = static_cast<PlayerState>(player_state);

    if (ImGui::Button("Step##SimPlayer", ImVec2(-1, 0))) {
        m_player_state = PlayerState::Playing;
        pre_draw_loop();
        m_player_state = PlayerState::Paused;
    }
    // --------------------------------------------------------------------
    ImGui::Checkbox("pause if intersecting", &m_bkp_has_intersections);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - stop playing if step has intersections.");

    ImGui::Checkbox("solve collisions", &m_state.m_solve_collisions);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - solve collisions automatically on each step.");

    ImGui::Checkbox("pause if optimization failed", &m_bkp_optimization_failed);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - stop playing if step optimization failed.");

    ImGui::Checkbox("pause on collisions", &m_bkp_had_collision);
    ImGui::SameLine();
    ImGui::HelpMarker("yes - stop playing if step had a collision.");

    // --------------------------------------------------------------------
    if (ImGui::SliderInt(
            "step##Replay", &m_state.m_num_simulation_steps, 0,
            m_state.state_sequence.size() - 1)) {
        replaying = true;
    }
}

void UISimState::draw_settings()
{
    auto config = m_state.get_active_config();
    ImGui::BeginChild(
        "##config", ImVec2(ImGui::GetWindowContentRegionWidth(), 300), false,
        ImGuiWindowFlags_HorizontalScrollbar);
    {
        ImGui::TreeNodeJson(config);
    }
    ImGui::EndChild();
}

void UISimState::draw_legends()
{
    static float second_col = 100;
    float slider_width = ImGui::GetWindowWidth() * 0.2f;
    for (auto& data : datas_) {
        const std::string& label = data.first;
        const auto& ptr = data.second;

        /// visibility checkbox
        bool tmp = ptr->visibility();
        if (!ptr->is_com()) {
            if (ImGui::Checkbox(("##UI-" + label).c_str(), &tmp)) {
                ptr->visibility(tmp);
            }
            /// Color
            ImGui::SameLine();
            if (ImGui::DoubleColorEdit3(
                    (label + "##UI-color").c_str(), ptr->m_color)) {
                ptr->recolor();
            }
        } else if (ImGui::Checkbox((label + "##UI-check").c_str(), &tmp)) {
            ptr->visibility(tmp);
        }

        /// other specific attributes
        if (ptr->is_mesh()) {
            ImGui::SameLine(second_col);
            if (ImGui::Checkbox(
                    ("data##UI-" + label).c_str(), &m_show_vertex_data)) {
            }
            ptr->show_vertex_data = m_show_vertex_data;
            ptr->update_vertex_data();

            ImGui::SameLine();
            ImGui::PushItemWidth(slider_width);
            {
                float point_size = ptr->data().point_size / pixel_ratio();
                if (ImGui::SliderFloat(
                        ("vertices##UI-scaling" + label).c_str(), &point_size,
                        0.00f, 10.0f, "%1.f")) {
                    ptr->data().point_size = point_size * pixel_ratio();
                }
            }
            ImGui::PopItemWidth();

        } else if (ptr->is_scalar_field()) {
            ImGui::SameLine();
            if (ImGui::Checkbox(
                    ("log##UI-" + label).c_str(), &ptr->m_use_log_scale)) {
                ptr->recolor();
            }
            ImGui::SameLine();
            bool show_iso = ptr->data().show_overlay;
            if (ImGui::Checkbox(("isolines##UI-" + label).c_str(), &show_iso)) {
                ptr->data().show_overlay = show_iso;
            }
            ImGui::SameLine();
            bool show_faces = ptr->data().show_faces;
            if (ImGui::Checkbox(("faces##UI-" + label).c_str(), &show_faces)) {
                ptr->data().show_faces = show_faces;
            }

        } else if (ptr->is_vector_field()) {
            ImGui::SameLine(second_col);
            if (ImGui::Checkbox(
                    ("norm##UI-" + label).c_str(), &ptr->m_normalized)) {
                ptr->recolor();
            }
            ImGui::SameLine();
            ImGui::PushItemWidth(slider_width);
            {
                float aux = float(ptr->m_scaling);
                if (ImGui::SliderFloat(
                        ("scale##UI-scaling" + label).c_str(), &aux, 0.00f,
                        10.0f, "%1.f")) {
                    ptr->m_scaling = double(aux);
                    ptr->recolor();
                }
            }
            ImGui::PopItemWidth();
        } else if (ptr->is_com()) {
            ImGui::SameLine();
            ImGui::PushItemWidth(ImGui::GetWindowWidth() - 150);
            {
                float aux = float(ptr->m_scaling);
                if (ImGui::SliderFloat(
                        ("scale##UI-com-scaling" + label).c_str(), &aux, 0.00f,
                        10.0f, "%1.1f", 1.0f)) {
                    ptr->m_scaling = double(aux);
                    ptr->recolor();
                }
            }
            ImGui::PopItemWidth();
        }
    }
}

} // namespace ipc::rigid
