#include <opt/solver_settings.hpp>

namespace ccd {
namespace opt {

    SolverSettings::SolverSettings(const OptimizationMethod method,
        const int verbosity, const int max_iter,
        const double relative_tolerance, const double absolute_tolerance,
        const double max_time, const callback_intermediate intermediate_cb,
        const QPSolver qp_solver, const NcpUpdate ncp_update_method,
        const LCPSolver lcp_solver, const double barrier_epsilon,
        const double min_barrier_epsilon, const double line_search_tolerance)
        : method(method)
        , verbosity(verbosity)
        , max_iter(max_iter)
        , relative_tolerance(relative_tolerance)
        , absolute_tolerance(absolute_tolerance)
        , max_time(max_time)
        , intermediate_cb(intermediate_cb)
        , qp_solver(qp_solver)
        , ncp_update_method(ncp_update_method)
        , lcp_solver(lcp_solver)
        , barrier_epsilon(barrier_epsilon)
        , min_barrier_epsilon(min_barrier_epsilon)
        , line_search_tolerance(line_search_tolerance)
    {
        if (this->intermediate_cb == nullptr) {
            this->intermediate_cb
                = [](const Eigen::VectorXd&, const double,
                      const Eigen::VectorXd&, const double, const int) {};
        }
    }

} // namespace opt
} // namespace ccd
