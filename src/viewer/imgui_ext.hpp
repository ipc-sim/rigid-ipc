#pragma once

#include <imgui/imgui.h>
#include <limits>

#define CCD_IM_ARRAYSIZE(_ARR) (int(sizeof(_ARR) / sizeof(*_ARR)))

namespace ImGui {

bool InputIntBounded(const char* label, int* val,
    int lower_bound = std::numeric_limits<int>::min(),
    int upper_bound = std::numeric_limits<int>::max(), int step = 1,
    int step_fast = 100, ImGuiInputTextFlags flags = 0);

bool InputDoubleBounded(const char* label, double* val,
    double lower_bound = -std::numeric_limits<double>::infinity(),
    double upper_bound = std::numeric_limits<double>::infinity(),
    double step = 1, double step_fast = 100, const char* format = "%.6f",
    ImGuiInputTextFlags flags = 0);

} // namespace ImGui
