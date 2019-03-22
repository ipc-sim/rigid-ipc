#include <opt/solver.hpp>

#include <opt/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/ipopt_solver.hpp>
#endif
#include <opt/linearized_constraint_solver.hpp>

#include <ccd/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    OptimizationResults solve_problem(
        OptimizationProblem& problem, const SolverSettings& settings)
    {
        if (!problem.validate_problem()) {
            return OptimizationResults();
        }

        switch (settings.method) {
        case MMA:
        case SLSQP:
            return solve_problem_with_nlopt(problem, settings);
        case IPOPT:
#ifdef BUILD_WITH_IPOPT
            return minimize_ipopt(problem, settings);
#else
            throw NotImplementedError("IPOPT is not enabled");
#endif
        case LINEARIZED_CONSTRAINTS:
            return solve_problem_with_linearized_constraints(problem, settings);
        case NCP:
            throw NotImplementedError(
                "Nonlinear complementarity problem not implemented");
        }
    }

} // namespace opt
} // namespace ccd
