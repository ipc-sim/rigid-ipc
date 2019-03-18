#pragma once

#include <Eigen/Core>

namespace ccd {
namespace opt {

    struct OptimizationResults {
        Eigen::MatrixXd x; ///< @brief the solution of the optimization.
        double minf;       ///< @brief value of the objective function
        bool success;      ///< @brief whether or not the optimizer exited
                           ///< successfully.

        OptimizationResults();
        OptimizationResults(Eigen::MatrixXd x, double minf, bool success);
    };

} // namespace opt
} // namespace ccd
