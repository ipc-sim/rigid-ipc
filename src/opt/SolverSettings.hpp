#pragma once

#include <Eigen/Core>
#include <opt/lcp_solver.hpp>
#include <opt/ncp_solver.hpp>

namespace ccd {
namespace opt {

    /**
     * @brief function type for callback, called once on each iteration
     * @param x             : current iterate of the optimization
     * @param obj_value     : value of the objective function
     * @param dual          : value of the dual variables
     * @param gamma         : last line-search step size
     * @param iteration     : current iteration number
     */
    typedef std::function<void(const Eigen::VectorXd& x, const double obj_value,
        const Eigen::VectorXd& dual, const double gamma, const int iteration)>
        callback_intermediate;

    /// @brief Methods available for non-linear constrained optimization.
    enum OptimizationMethod {
        MMA,   ///< @brief Method of Moving Asymptotes (NLopt)
        SLSQP, ///< @brief Sequential Least-Squares Quadratic Programming
               ///< (NLopt)
        IPOPT, ///< @brief Interior-Point Method (Ipopt)
        LINEARIZED_CONSTRAINTS, ///< @brief Linearize the constraints and solve
                                ///< the QP (OSQP/MOSEK)
        NCP,                    ///< @brief Nonlinear Complementarity Problem
        BARRIER_NEWTON          ///< @brief Barrier Newton's Method
    };

    static const char* OptimizationMethodNames[] = { "MMA", "SLSQP", "IPOPT",
        "Linearized Const.", "NCP", "Barrier Newton" };

    enum QPSolver {
        OSQP, ///< @brief Use OSQP to solve the qudratic program
        MOSEK ///< @brief Use MOSEK to solve the qudratic program
    };

    static const char* QPSolverNames[] = { "OSQP", "MOSEK" };

    struct SolverSettings {
        OptimizationMethod method;   ///< @brief Optimization method to use
        QPSolver qp_solver;          ///< @brief Quadratic programming solver
        LCPSolver lcp_solver;        ///< @brief LCP solver
        NcpUpdate ncp_update_method; ///< @brief Method to update NCP solution
        int verbosity;               ///< @brief Verbosity of output
        int max_iter;                ///< @brief Maximum number of iterations
        double relative_tolerance;   ///< @brief Relative tolerance for x
        double absolute_tolerance;   ///< @brief Absolute tolerance for x
        double constraint_tolerance; ///< @brief Tolerance of the constraint
        double max_time;             ///< @brief Max time to spend solving
        /// @brief Callback function, called every outer iteration
        callback_intermediate intermediate_cb;

        // ToDo: Find a better place to store this
        /// @brief Starting epsilon for barrier methods using the spline_barrier
        double starting_barrier_epsilon;

        /// @brief Construct the solver settings
        SolverSettings(const OptimizationMethod method = SLSQP,
            const QPSolver qp_solver = OSQP, const int verbosity = 0,
            const int max_iter = 3000, const double relative_tolerance = 1e-8,
            const double absolute_tolerance = 1e-8,
            const double constraint_tolerance = 1e-8,
            const double max_time = 2e19,
            const callback_intermediate intermediate_cb = nullptr,
            const double starting_barrier_epsilon = 0.0);
    };

} // namespace opt
} // namespace ccd
