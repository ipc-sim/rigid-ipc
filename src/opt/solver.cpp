#include <opt/solver.hpp>

#include <opt/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/ipopt_solver.hpp>
#endif
#include <opt/barrier_newton_solver.hpp>
#include <opt/linearized_constraint_solver.hpp>

#include <ccd/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    OptimizationResults solve_problem(
        OptimizationProblem& problem, SolverSettings& settings)
    {
        if (!problem.validate_problem()) {
            return OptimizationResults();
        }

        switch (settings.method) {
        case OptimizationMethod::BARRIER_NEWTON:
            return solve_problem_with_barrier_newton(problem, settings);
        default:
            throw NotImplementedError("Not implemented");
        }
    }

    OptimizationSolver::OptimizationSolver() {}
    OptimizationSolver::~OptimizationSolver() {}

} // namespace opt
} // namespace ccd
