#pragma once

#include <Eigen/Core>

namespace ccd {
namespace opt {

    ///@brief function type for callback, called once on each iteration
    typedef std::function<void(const Eigen::VectorXd& x, const double obj_value,
        const Eigen::VectorXd& dual, const int iteration)>
        callback_intermediate;

    /// @brief Methods available for non-linear constrained optimization.
    enum OptimizationMethod {
        MMA,   ///< @brief Method of Moving Asymptotes (NLopt)
        SLSQP, ///< @brief Sequential Least-Squares Quadratic Programming
               ///< (NLopt)
        IP,    ///< @brief Interior-Point Method (Ipopt)
        LINEARIZED_CONSTRAINTS, ///< @brief Linearize the constraints and solve
                                ///< the QP (OSQP/Mosek)
        NCP                     ///< @brief Nonlinear Complementarity Problem
    };

    static const char* OptimizationMethodStrings[]
        = { "mma", "slsqp", "ipopt", "linearized constraints", "ncp" };

    struct SolverSettings {
        OptimizationMethod method;   ///< @brief Optimization method to use
        int verbosity;               ///< @brief Verbosity of output
        int max_iter;                ///< @brief Maximum number of iterations
        double relative_tolerance;   ///< @brief Relative tolerance for x
        double absolute_tolerance;   ///< @brief Absolute tolerance for x
        double constraint_tolerance; ///< @brief Tolerance of the constraint
        double max_time;             ///< @brief Max time to spend solving
        callback_intermediate intermediate_cb; ///< @brief Callback function,
                                               ///< called every outer iteration

        /// @brief Default constructor
        SolverSettings();

        /// @brief Construct the solver settings
        SolverSettings(const OptimizationMethod method, const int verbosity = 0,
            const int max_iter = 3000, const double relative_tolerance = 1e-8,
            const double absolute_tolerance = 1e-8,
            const double constraint_tolerance = 1e-8,
            const double max_time = 2e19,
            const callback_intermediate intermediate_cb = nullptr);
    };

} // namespace opt
} // namespace ccd
