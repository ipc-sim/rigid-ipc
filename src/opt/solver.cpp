#include <opt/solver.hpp>

#include <opt/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/ipopt_solver.hpp>
#endif
#ifdef BUILD_WITH_OSQP
#include <opt/linearized_constraint_solver.hpp>
#endif
// #include <opt/ncp_solver.hpp> // TODO: Create this file

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
            // Implemented in NLopt
            return solve_problem_with_nlopt(problem, settings);
        case IP:
#ifdef BUILD_WITH_IPOPT
            // Implemented in Ipopt
            return ccd::opt::minimize_ipopt(problem);
#else
            throw NotImplementedError("IPOPT not Enabled");
#endif
        case LINEARIZED_CONSTRAINTS:
#ifdef BUILD_WITH_OSQP
            return solve_problem_with_linearized_constraints(problem, settings);
#else
            throw NotImplementedError("OSQP not Enabled");
#endif
        case NCP:
            throw NotImplementedError(
                "Nonlinear complementarity problem not implemented");
        }
    }

} // namespace opt
} // namespace ccd
