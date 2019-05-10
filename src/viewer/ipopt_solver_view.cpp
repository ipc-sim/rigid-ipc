#include "ipopt_solver_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void ipopt_solver_menu(ccd::opt::IpoptSolver& solver)
{
    ImGui::InputIntBounded("max iter##ipopt_solver", &solver.max_iterations, 0);
    ImGui::InputIntBounded("verbose##ipopt_solver", &solver.print_level, 0);
    ImGui::InputDoubleBounded("tolerance##ipopt_solver", &solver.tolerance, 0.0,
        2e19, 0.0, 0.0, "%.3g");
}
} // namespace ccd
