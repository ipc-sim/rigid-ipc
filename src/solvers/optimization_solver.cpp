#include "optimization_solver.hpp"

#include <solvers/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <solvers/ipopt_solver.hpp>
#endif
#include <solvers/barrier_newton_solver.hpp>
#include <solvers/qp_solver.hpp>

#include <ccd/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    OptimizationSolver::OptimizationSolver() {}
    OptimizationSolver::~OptimizationSolver() {}

} // namespace opt
} // namespace ccd
