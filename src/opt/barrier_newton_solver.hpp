#pragma once

/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <opt/optimization_solver.hpp>
#include <opt/barrier_constraint.hpp>

namespace ccd {
namespace opt {


    /**
     * @brief Apply the constraints of the general problem as a barrier added to
     * the objective of the original problem.
     *
     * @param[in] general_problem Problem to convert to a problem with barrier
     *                            on the objective
     * @param[in] epsilon The barrier epsilon controlling the steepness. This
     *                    value is stored by reference in the barrier functions.
     * @paramp[out] barrier_problem Problem to store the modified objective of
     *                              the general problem.
     */

    class BarrierProblem : public OptimizationProblem {
    public:
        BarrierProblem(OptimizationProblem& problem, double epsilon);
        ~BarrierProblem() override;

        Eigen::VectorXd barrier(const Eigen::VectorXd x);
        Eigen::VectorXd barrier_gradient(const Eigen::VectorXd x);
        Eigen::VectorXd barrier_hessian(const Eigen::VectorXd x);

        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_g(const Eigen::VectorXd&) override
        {
            return Eigen::VectorXd();
        }
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd&) override
        {
            return Eigen::MatrixXd();
        }
        std::vector<Eigen::MatrixXd> eval_hessian_g(
            const Eigen::VectorXd&) override
        {
            return std::vector<Eigen::MatrixXd>();
        }
        OptimizationProblem* general_problem;
        double epsilon;
    };

    class BarrierNewtonSolver : public OptimizationSolver {
    public:
        BarrierNewtonSolver();
        ~BarrierNewtonSolver() override;
        OptimizationResults solve(OptimizationProblem& problem) override;

        BarrierConstraint * barrier_constraint;

        double min_barrier_epsilon;
        double absolute_tolerance;
        double line_search_tolerance;
        int max_iterations;
    };

} // namespace opt
} // namespace ccd
