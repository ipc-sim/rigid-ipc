#include <opt/SolverSettings.hpp>

namespace ccd {
namespace opt {

    SolverSettings::SolverSettings(const OptimizationMethod method,
        const QPSolver qp_solver, const int verbosity, const int max_iter,
        const double relative_tolerance, const double absolute_tolerance,
        const double constraint_tolerance, const double max_time,
        const callback_intermediate intermediate_cb)
        : method(method)
        , qp_solver(qp_solver)
        , lcp_solver(LCP_GAUSS_SEIDEL)
        , ncp_update_method(NcpUpdate::G_GRADIENT)
        , verbosity(verbosity)
        , max_iter(max_iter)
        , relative_tolerance(relative_tolerance)
        , absolute_tolerance(absolute_tolerance)
        , constraint_tolerance(constraint_tolerance)
        , max_time(max_time)
        , intermediate_cb(intermediate_cb)
    {
    }

} // namespace opt
} // namespace ccd
