#include "imgui_ext.hpp"

#include <imgui.h>

namespace ImGui {

void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void TreeNodeJson(const nlohmann::json json)
{
    static const ImVec4 label_color =
        ImVec4(241.f / 255.f, 196.f / 255.f, 15.f / 255.f, 1.0f);
    ImGui::Unindent(ImGui::GetTreeNodeToLabelSpacing() * 0.1f);
    ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 1);
    for (auto& j : json.items()) {
        if (j.value().is_object()) {
            ImGui::PushStyleColor(ImGuiCol_Text, label_color);
            if (ImGui::TreeNode(j.key().c_str())) {
                ImGui::PopStyleColor();
                TreeNodeJson(j.value());
                ImGui::TreePop();
            } else {
                ImGui::PopStyleColor();
            }
        } else {
            ImGui::Bullet();
            ImGui::TextColored(label_color, "%s: ", j.key().c_str());
            ImGui::SameLine();
            if (j.value().is_number_float()) {
                ImGui::Text("%g", j.value().get<double>());
            } else {
                ImGui::Text("%s", j.value().dump().c_str());
            }
        }
    }
    ImGui::PopStyleVar();
    ImGui::Indent(ImGui::GetTreeNodeToLabelSpacing() * 0.1f);
}

} // namespace ImGui
