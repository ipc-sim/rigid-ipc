#include "constraint_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void collision_constraint_menu(ccd::opt::CollisionConstraint& cstr)
{
    ImGui::Checkbox("update col. set##opt", &(cstr.update_collision_set));
    ImGui::Checkbox("extend col. set##opt", &(cstr.extend_collision_set));

    int idx_detection_method = cstr.detection_method;
    if (ImGui::Combo("method##ccd", &idx_detection_method,
            ccd::DetectionMethodNames,
            CCD_IM_ARRAYSIZE(ccd::DetectionMethodNames))) {
        cstr.detection_method
            = static_cast<ccd::DetectionMethod>(idx_detection_method);
    }
}

void volume_constraint_menu(ccd::opt::VolumeConstraint& cstr)
{
    collision_constraint_menu(cstr);
    ImGui::InputDoubleBounded(
        "volume eps.##opt", &cstr.volume_epsilon, 0.0, 2e19, 0.0, 0.0, "%.3g");
}

void barrier_constraint_menu(ccd::opt::BarrierConstraint& cstr)
{
    collision_constraint_menu(cstr);
    int initial_epsilon_idx = static_cast<int>(cstr.initial_epsilon);
    if (ImGui::Combo("init. eps##opt", &initial_epsilon_idx,
            ccd::opt::InitialBarrierEpsilonNames,
            CCD_IM_ARRAYSIZE(ccd::opt::InitialBarrierEpsilonNames))) {
        cstr.initial_epsilon
            = static_cast<ccd::opt::InitialBarrierEpsilon>(initial_epsilon_idx);
    }
}

} // namespace ccd
