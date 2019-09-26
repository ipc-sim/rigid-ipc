#pragma once

#include <Eigen/Core>

namespace ccd {
namespace opt {

    struct OptimizationResults {
        Eigen::MatrixXd x; ///< @brief The solution of the optimization.
        double minf;       ///< @brief Value of the objective function.
        bool success;      ///< @brief Whether or not the optimizer exited
                           ///< successfully.
        bool finished;     ///< @brief Whether or not the optimizer converged.

        OptimizationResults();
        OptimizationResults(
            Eigen::MatrixXd x, double minf, bool success, bool finished);
    };

} // namespace opt
} // namespace ccd
