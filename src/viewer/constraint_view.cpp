#include "constraint_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void collision_constraint_menu(ccd::opt::CollisionConstraint& cstr)
{
    ImGui::Checkbox("recompute col. set##opt", &(cstr.recompute_collision_set));
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
}

} // namespace ccd
