#pragma once

#include <nlohmann/json.hpp>

#define CCD_IM_ARRAYSIZE(_ARR) (int(sizeof(_ARR) / sizeof(*_ARR)))

namespace ImGui {

void TreeNodeJson(const nlohmann::json json);
void HelpMarker(const char* desc);

} // namespace ImGui
