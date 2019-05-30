#include "imgui_ext.hpp"

namespace ImGui {

bool InputIntBounded(const char* label, int* val, int lower_bound,
    int upper_bound, int step, int step_fast, ImGuiInputTextFlags flags)
{
    int unbounded_val = *val;
    if (ImGui::InputInt(label, &unbounded_val, step, step_fast, flags)) {
        if (unbounded_val >= lower_bound && unbounded_val <= upper_bound) {
            *val = unbounded_val;
            return true;
        }
    }
    return false;
};

bool InputDoubleBounded(const char* label, double* val, double lower_bound,
    double upper_bound, double step, double step_fast, const char* format,
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
};

} // namespace ImGui
