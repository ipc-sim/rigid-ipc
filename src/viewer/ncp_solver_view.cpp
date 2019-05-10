#include <viewer/ncp_solver_view.hpp>

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void ncp_solver_menu(ccd::opt::NCPDisplacementOptimization& ncp_solver)
{
    int idx_lcp_solver = ncp_solver.lcp_solver;
    int idx_ncp_update = static_cast<int>(ncp_solver.update_method);

    ImGui::InputIntBounded(
        "max iter##ncp_solver", &ncp_solver.max_iterations, 0);

    if (ImGui::Combo("x-update##opt", &idx_ncp_update,
            ccd::opt::NcpUpdateNames,
            CCD_IM_ARRAYSIZE(ccd::opt::NcpUpdateNames))) {
        ncp_solver.update_method
            = static_cast<ccd::opt::NcpUpdate>(idx_ncp_update);
    }

    if (ImGui::Combo("LCP solver##opt", &idx_lcp_solver,
            ccd::opt::LCPSolverNames,
            CCD_IM_ARRAYSIZE(ccd::opt::LCPSolverNames))) {
        ncp_solver.lcp_solver
            = static_cast<ccd::opt::LCPSolver>(idx_lcp_solver);
    }

    ImGui::Checkbox("keep unf.##ncp_solver", &ncp_solver.keep_in_unfeasible);
    ImGui::Checkbox("##ncp_solver", &ncp_solver.check_convergence);
    ImGui::SameLine();
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
    ImGui::InputDoubleBounded("conv. tol.##ncp_solver",
        &ncp_solver.convegence_tolerance, 0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::PopItemWidth();
}
} // namespace ccd
