#pragma once

/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */
#include <memory>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

    class BarrierProblem : public OptimizationProblem {
    public:
        BarrierProblem(OptimizationProblem& problem);
        ~BarrierProblem() override {}

        double eval_f(const Eigen::VectorXd& x) override;

        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;

        /////
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

        /////
        void enable_line_search_mode(const Eigen::VectorXd& max_x) override;
        void disable_line_search_mode() override;

        bool eval_intermediate_callback(const Eigen::VectorXd& x) override;

        bool has_barrier_constraint() override
        {
            return general_problem->has_barrier_constraint();
        }
        double get_barrier_epsilon() override
        {
            return general_problem->get_barrier_epsilon();
        }

        const Eigen::VectorXb& is_dof_fixed() override;
        OptimizationProblem* general_problem;
    };

    class BarrierSolver : public OptimizationSolver {
    public:
        BarrierSolver();
        BarrierSolver(const std::string& name);
        ~BarrierSolver() override {}

        OptimizationResults solve(OptimizationProblem& problem) override;
        OptimizationSolver& get_inner_solver() const
        {
            return *inner_solver_ptr;
        }

        void clear() override;

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        void init(OptimizationProblem& problem) override;
        OptimizationResults step_solve() override;
        int num_outer_iterations() override { return num_outer_iterations_; }
        inline double barrier_epsilon()
        {
            return general_problem_ptr->get_barrier_epsilon();
        }

        void eval_f(
            const Eigen::MatrixXd& points, Eigen::VectorXd& fx) override;

        Eigen::VectorXd get_grad_f() const override;
        bool has_inner_solver() override { return true; }
        const OptimizationSolver& inner_solver() override
        {
            return get_inner_solver();
        }

        double min_barrier_epsilon;

    protected:
        std::shared_ptr<OptimizationSolver> inner_solver_ptr;

        std::unique_ptr<BarrierProblem> barrier_problem_ptr;
        OptimizationProblem* general_problem_ptr;
        int num_outer_iterations_;
    };

} // namespace opt
} // namespace ccd
