#include "viewer.hpp"

namespace ccd {
void ViewerMenu::draw_menu()
{
    float menu_width = 180.f * menu_scaling();
    draw_labels_window();

    ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    bool _viewer_menu_visible = true;
    ImGui::Begin(
        "CCD Viewer", &_viewer_menu_visible, ImGuiWindowFlags_NoSavedSettings  | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);

    draw_io();
    draw_edit_modes();
    draw_ui_settings();
    draw_ccd_steps();

    ImGui::PopItemWidth();
    ImGui::End();


    ImGui::SetNextWindowPos(ImVec2(menu_width + 10.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(800.0f, 50.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(200.0f, 30.0f), ImVec2(800.0f, -1.0f));
    bool _message_menu_visible = true;
    ImGui::Begin(
        "CCD Message", &_message_menu_visible, ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    static ImVec4 color_ok(0.0,1.0,0.0,1.0);
    static ImVec4 color_error(1.0,0.0,0.0,1.0);

    ImGui::TextColored(last_action_success ? color_ok: color_error, "%s", last_action_message.c_str());
    ImGui::End();
}

void ViewerMenu::draw_io(){

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

void color_button(){
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.5,0.0,0.0,1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f,0.0,0.0,1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(1.0,0.0,0.0,1.0));
}

void ViewerMenu::draw_edit_modes(){

    if (ImGui::CollapsingHeader("Edit Mode", ImGuiTreeNodeFlags_DefaultOpen)) {
        float w = ImGui::GetContentRegionAvailWidth();
        float p = ImGui::GetStyle().FramePadding.x;

        for(uint i=0; i < ViewerEditModeAll.size(); ++i){
            bool needs_pop = false;
            if (edit_mode == ViewerEditModeAll[i]){
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
            if (i%2 == 0){
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
    if (ImGui::CollapsingHeader("UI Settings", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Checkbox("show surface", &(viewer->data_list[surface_data_id].show_overlay));
        ImGui::Checkbox("show vertex id", &(viewer->data_list[surface_data_id].show_vertid));
        ImGui::InputFloat("point size ##surface", &(viewer->data_list[surface_data_id].point_size));

        ImGui::Checkbox("show displacements", &(viewer->data_list[displ_data_id].show_overlay));
        ImGui::InputFloat("point size ##displ", &(viewer->data_list[displ_data_id].point_size));
        ImGui::Checkbox("show gradient", &(viewer->data_list[gradient_data_id].show_overlay));

        if(ImGui::SliderFloat("time", &(state.time), 0.0, 1.0)){
            redraw_at_time();
        }
    }
}

void ViewerMenu::draw_ccd_steps()
{
    static int idx_detection_method = 0;
    if (ImGui::CollapsingHeader("CCD Steps", ImGuiTreeNodeFlags_DefaultOpen)) {
        std::string combo_options;
        for (const auto &piece : DetectionMethodNames) combo_options += piece + "\0";
        combo_options+="\0";

        if (ImGui::Combo("Method", &idx_detection_method, combo_options.c_str())){
            state.detection_method = DetectionMethodAll[size_t(idx_detection_method)];
        }
        if (ImGui::Button("detect EV collisions", ImVec2(-1,0))){
            detect_edge_vertex_collisions();
        }

        if(state.impacts != nullptr && state.impacts->size()){
            int impact_id = state.current_impact + 1;
            ImGui::Text("impact %i/%i", state.current_impact+1, int(state.impacts->size()));

            if(ImGui::InputInt("", &impact_id)){
                goto_impact(impact_id - 1);
            }

            if (ImGui::Button("Goto Next impact", ImVec2(-1,0))){
                goto_impact(state.current_impact + 1);
            }
        }
        ImGui::Separator();
        ImGui::InputDouble("vol. epsilon", &state.epsilon);
    }

}
}
