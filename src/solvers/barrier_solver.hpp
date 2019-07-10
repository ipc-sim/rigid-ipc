#pragma once

/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */

#include <opt/barrier_constraint.hpp>
#include <solvers/bfgs_solver.hpp>
#include <solvers/gradient_descent_solver.hpp>
#include <solvers/newton_solver.hpp>

namespace ccd {
namespace opt {

    enum class BarrierInnerSolver { NEWTON, BFGS, GRADIENT_DESCENT };
    static const char* BarrierInnerSolverNames[]
        = { "Newton's Method", "BFGS", "Gradient Descent" };

    class BarrierSolver : public OptimizationSolver {
    public:
        BarrierSolver();
        ~BarrierSolver() override;
        OptimizationResults solve(OptimizationProblem& problem) override;
        OptimizationSolver& get_inner_solver();

        double min_barrier_epsilon;
        BarrierInnerSolver inner_solver_type;

    protected:
        NewtonSolver newton_inner_solver;
        BFGSSolver bfgs_inner_solver;
        GradientDescentSolver gradient_descent_inner_solver;
    };

    class BarrierProblem : public OptimizationProblem {
    public:
        BarrierProblem(OptimizationProblem& problem, double epsilon);
        ~BarrierProblem() override;

        Eigen::VectorXd barrier(const Eigen::VectorXd x);
        Eigen::VectorXd barrier_gradient(const Eigen::VectorXd x);
        Eigen::VectorXd barrier_hessian(const Eigen::VectorXd x);

        double eval_f(const Eigen::VectorXd& x) override;

        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x, double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;

        Eigen::VectorXd eval_g(const Eigen::VectorXd&) override
        {
            return Eigen::VectorXd();
        }

        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd&) override
        {
            return Eigen::MatrixXd();
        }

        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd&) override
        {
            return std::vector<Eigen::SparseMatrix<double>>();
        }

        void enable_line_search_mode(const Eigen::VectorXd& max_x) override;
        void disable_line_search_mode() override;

        bool eval_intermediate_callback(const Eigen::VectorXd& x) override;
        OptimizationProblem* general_problem;

        double epsilon;
    };

} // namespace opt
} // namespace ccd
