#pragma once

#include <Eigen/Core>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {
    enum class QPSolverType {
        OSQP, ///< @brief Use OSQP to solve the qudratic program
        MOSEK ///< @brief Use MOSEK to solve the qudratic program
    };

    class QPSolver : public IOptimizationSolver {
    public:
        QPSolver();
        ~QPSolver() override;

        OptimizationResults solve(OptimizationProblem& problem) override;

        double absolute_tolerance;
        double relative_tolerance;
        int max_iterations;

        QPSolverType qp_solver;

        bool solve_with_osqp(Eigen::MatrixXd& x);

        bool solve_with_mosek(Eigen::MatrixXd& x);

        // Optimization Problem structure
        // Minimize: ½ * x^T * Q * x + c^T * x + cf
        //
        // Subject to: lc ≤ Ax ≤ uc
        //             lx ≤ x ≤ ux
        Eigen::SparseMatrix<double> Q;
        Eigen::VectorXd c;
        double cf;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd lc, uc;
        Eigen::VectorXd lx, ux;
        void quadratic_energy(OptimizationProblem& problem);
        void linearized_constraints(OptimizationProblem& problem);
    };

} // namespace opt

} // namespace ccd
