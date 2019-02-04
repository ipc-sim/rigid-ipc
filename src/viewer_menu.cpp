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

    ImGui::PopItemWidth();
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

        if(ImGui::SliderFloat("displ", &(state.displacement_ptge), 0.0, 1.0)){
            redraw_with_displacements();
        }
    }
}
}
