#pragma once

#include <Eigen/Core>

namespace ipc::rigid {

struct OptimizationResults {
    Eigen::VectorXd x;  ///< @brief The solution of the optimization.
    double minf;        ///< @brief Value of the objective function.
    bool success;       ///< @brief Whether or not the optimizer exited
                        ///< successfully.
    bool finished;      ///< @brief Whether or not the optimizer converged.
    int num_iterations; ///< @brief number of iterations performed

    OptimizationResults();
    OptimizationResults(
        const Eigen::VectorXd& x,
        double minf,
        bool success,
        bool finished,
        int num_iterations);
};

} // namespace ipc::rigid
