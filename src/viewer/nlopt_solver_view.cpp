#include "nlopt_solver_view.hpp"

#include <viewer/ncp_solver_view.hpp>

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {
void nlopt_solver_menu(ccd::opt::NLOptSolver& solver)
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
