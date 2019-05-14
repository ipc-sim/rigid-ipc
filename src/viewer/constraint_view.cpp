#include "constraint_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void volume_constraint_menu(ccd::opt::VolumeConstraint& cstr)
{
    ImGui::InputDoubleBounded(
        "volume eps.##opt", &cstr.volume_epsilon, 0.0, 2e19, 0.0, 0.0, "%.3g");
}

void barrier_constraint_menu(ccd::opt::BarrierConstraint& /*cstr*/)
{

}

} // namespace ccd
