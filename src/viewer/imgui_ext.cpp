#include "imgui_ext.hpp"

namespace ImGui {

bool InputIntBounded(const char* label,
    int* val,
    int lower_bound,
    int upper_bound,
    int step,
    int step_fast,
    ImGuiInputTextFlags flags)
{
    int unbounded_val = *val;
    if (ImGui::InputInt(label, &unbounded_val, step, step_fast, flags)) {
        if (unbounded_val >= lower_bound && unbounded_val <= upper_bound) {
            *val = unbounded_val;
            return true;
        }
    }
    return false;
}

bool InputDoubleBounded(const char* label,
    double* val,
    double lower_bound,
    double upper_bound,
    double step,
    double step_fast,
    const char* format,
    ImGuiInputTextFlags flags)
{
    double unbounded_val = *val;
    if (ImGui::InputDouble(
            label, &unbounded_val, step, step_fast, format, flags)) {
        if (unbounded_val >= lower_bound && unbounded_val <= upper_bound) {
            *val = unbounded_val;
            return true;
        }
    }
    return false;
}

bool DragDouble(const char* label,
    double* v,
    double v_speed,
    double v_min,
    double v_max,
    const char* format,
    float power)
{

    return DragScalar(
        label, ImGuiDataType_Double, v, float(v_speed), &v_min, &v_max, format, power);
}

bool DoubleColorEdit3(const char* label, Eigen::RowVector3d& color)
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

void HelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

} // namespace ImGui
