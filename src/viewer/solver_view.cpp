#include "solver_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {

void solver_menu(ccd::opt::QPSolver& solver)
{
    int idx_algorithm = static_cast<int>(solver.qp_solver);

    ImGui::InputIntBounded("max iter##qp_solver", &solver.max_iterations, 0);

    if (ImGui::Combo("algorithm##qp_solver", &idx_algorithm,
            ccd::opt::QPSolverNames,
            CCD_IM_ARRAYSIZE(ccd::opt::QPSolverNames))) {
        solver.qp_solver = static_cast<ccd::opt::QPSolverType>(idx_algorithm);
    }

    ImGui::InputDoubleBounded("rel. tol.##nlopt", &solver.relative_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("abs. tol.##nlopt", &solver.absolute_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");
}

void solver_menu(ccd::opt::BarrierNewtonSolver& solver)
{
    ImGui::InputIntBounded(
        "max iter##barrier_solver", &solver.max_iterations, 0);
    ImGui::InputDoubleBounded("eps##barrier_solver", &solver.barrier_epsilon,
        0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("min eps##barrier_solver",
        &solver.min_barrier_epsilon, 0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("tol. abs##barrier_solver",
        &solver.absolute_tolerance, 0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("tol. line_search##barrier_solver",
        &solver.line_search_tolerance, 0.0, 2e19, 0.0, 0.0, "%.3g");
}

void solver_menu(ccd::opt::IpoptSolver& solver)
{
    ImGui::InputIntBounded("max iter##ipopt_solver", &solver.max_iterations, 0);
    ImGui::InputIntBounded("verbose##ipopt_solver", &solver.print_level, 0);
    ImGui::InputDoubleBounded("tolerance##ipopt_solver", &solver.tolerance, 0.0,
        2e19, 0.0, 0.0, "%.3g");
}

void solver_menu(ccd::opt::NCPSolver& ncp_solver)
{
    int idx_lcp_solver = ncp_solver.lcp_solver;
    int idx_ncp_update = static_cast<int>(ncp_solver.update_type);

    ImGui::InputIntBounded(
        "max iter##ncp_solver", &ncp_solver.max_iterations, 0);

    if (ImGui::Combo("x-update##opt", &idx_ncp_update, ccd::opt::NcpUpdateNames,
            CCD_IM_ARRAYSIZE(ccd::opt::NcpUpdateNames))) {
        ncp_solver.update_type
            = static_cast<ccd::opt::NcpUpdate>(idx_ncp_update);
    }

    if (ImGui::Combo("LCP solver##opt", &idx_lcp_solver,
            ccd::opt::LCPSolverNames,
            CCD_IM_ARRAYSIZE(ccd::opt::LCPSolverNames))) {
        ncp_solver.lcp_solver
            = static_cast<ccd::opt::LCPSolver>(idx_lcp_solver);
    }

    ImGui::Checkbox(
        "keep unfeasible##ncp_solver", &ncp_solver.keep_in_unfeasible);
    ImGui::Checkbox("##ncp_solver", &ncp_solver.check_convergence);
    ImGui::SameLine();
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
    ImGui::InputDoubleBounded("conv. tol.##ncp_solver",
        &ncp_solver.convergence_tolerance, 0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::PopItemWidth();
}

void solver_menu(ccd::opt::NLOptSolver& solver)
{
    int idx_algorithm;
    for (size_t i = 0; i < ccd::opt::NLOptAlgorithm.size(); i++) {
        if (solver.algorithm == ccd::opt::NLOptAlgorithm[i]) {
            idx_algorithm = int(i);
            break;
        }
    }

    ImGui::InputIntBounded("max iter##nlopt_solver", &solver.max_iterations, 0);

    if (ImGui::Combo("algorithm##nlopt_solver", &idx_algorithm,
            ccd::opt::NLOptAlgorithmNames,
            CCD_IM_ARRAYSIZE(ccd::opt::NLOptAlgorithmNames))) {
        solver.algorithm = ccd::opt::NLOptAlgorithm[size_t(idx_algorithm)];
    }

    ImGui::InputDoubleBounded("rel. tol.##nlopt", &solver.relative_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("abs. tol.##nlopt", &solver.absolute_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded(
        "max_time##nlopt", &solver.max_time, 0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::Checkbox("verbosity##opt", &solver.verbosity);
}
} // namespace ccd
