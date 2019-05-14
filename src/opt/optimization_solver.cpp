#include <opt/optimization_solver.hpp>

#include <opt/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/ipopt_solver.hpp>
#endif
#include <opt/barrier_newton_solver.hpp>
#include <opt/qp_solver.hpp>

#include <ccd/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    OptimizationSolver::OptimizationSolver() {}
    OptimizationSolver::~OptimizationSolver() {}

} // namespace opt
} // namespace ccd
