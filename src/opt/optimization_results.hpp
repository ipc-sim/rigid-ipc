#pragma once

#include <Eigen/Core>

#include <opt/solver_settings.hpp>

namespace ccd {
namespace opt {

    struct OptimizationResults {
        Eigen::MatrixXd x; ///< @brief the solution of the optimization.
        double minf;       ///< @brief value of the objective function
        bool success;      ///< @brief whether or not the optimizer exited
                           ///< successfully.
        OptimizationMethod method; ///<@brief method used to generate results
        bool finished = false;

        OptimizationResults();
        OptimizationResults(Eigen::MatrixXd x, double minf, bool success);
    };

} // namespace opt
} // namespace ccd
