#pragma once

#include <imgui/imgui.h>
#include <limits>

#define CCD_IM_ARRAYSIZE(_ARR) (int(sizeof(_ARR) / sizeof(*_ARR)))

namespace ImGui {

bool InputIntBounded(const char* label,
    int* val,
    int lower_bound = std::numeric_limits<int>::min(),
    int upper_bound = std::numeric_limits<int>::max(),
    int step = 1,
    int step_fast = 100,
    ImGuiInputTextFlags flags = 0);

bool InputDoubleBounded(const char* label,
    double* val,
    double lower_bound = -std::numeric_limits<double>::infinity(),
    double upper_bound = std::numeric_limits<double>::infinity(),
    double step = 1,
    double step_fast = 100,
    const char* format = "%.6f",
    ImGuiInputTextFlags flags = 0);

bool DragDouble(const char* label,
    double* v,
    double v_speed = 1.0,
    double v_min = 0.0,
    double v_max = 0.0,
    const char* format = "%.3f",
    float power = 1.0f); // If v_min >= v_max we have no bound

} // namespace ImGui
