#include "linearized_constraint_solver_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {

void linearized_constraint_solver_view(
    ccd::opt::LinearizedCstrSolver& solver)
{
    int idx_algorithm = static_cast<int>(solver.qp_solver);

    ImGui::InputIntBounded("max iter##qp_solver", &solver.max_iterations, 0);

    if (ImGui::Combo("algorithm##qp_solver", &idx_algorithm,
            ccd::opt::QPSolverNames,
            CCD_IM_ARRAYSIZE(ccd::opt::QPSolverNames))) {
        solver.qp_solver = static_cast<ccd::opt::QPSolver>(idx_algorithm);
    }

    ImGui::InputDoubleBounded("rel. tol.##nlopt", &solver.relative_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");
    ImGui::InputDoubleBounded("abs. tol.##nlopt", &solver.absolute_tolerance,
        0.0, 2e19, 0.0, 0.0, "%.3g");

}

} // namespace ccd
