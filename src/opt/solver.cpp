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
        case OptimizationMethod::MMA:
        case OptimizationMethod::SLSQP:
            return solve_problem_with_nlopt(problem, settings);
        case OptimizationMethod::IPOPT:
#ifdef BUILD_WITH_IPOPT
            return minimize_ipopt(problem, settings);
#else
            throw NotImplementedError("IPOPT is not enabled");
#endif
        case OptimizationMethod::LINEARIZED_CONSTRAINTS:
            return solve_problem_with_linearized_constraints(problem, settings);
        case OptimizationMethod::NCP:
            // ToDo: Call NCP from here instead of displacment_opt
            throw NotImplementedError(
                "Nonlinear complementarity problem not implemented");
        case OptimizationMethod::BARRIER_NEWTON:
            return solve_problem_with_barrier_newton(problem, settings);
        }
    }

} // namespace opt
} // namespace ccd
