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
        double eval_f(const Eigen::VectorXd& x,
            const bool update_constraint_set) override;

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

    class BarrierSolver : public IFullOptimizationSolver {
    public:
        BarrierSolver();
        BarrierSolver(const std::string& name);
        ~BarrierSolver() override {}

        // From IOptimizationSolver
        OptimizationResults solve(OptimizationProblem& problem) override;

        // From IFullOptimizationSolver
        void init(OptimizationProblem& problem) override;
        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        const std::string& name() const override { return name_; }
        bool has_inner_solver() override { return true; }
        const IBarrierOptimizationSolver& inner_solver() override
        {
            return *inner_solver_ptr;
        }

        OptimizationResults step_solve() override;
        Eigen::VectorXd get_grad_kkt() const override;
        int num_outer_iterations() const override
        {
            return num_outer_iterations_;
        }

        double min_barrier_epsilon;
        int max_iterations;

    protected:
        IBarrierOptimizationSolver& get_inner_solver() const
        {
            return *inner_solver_ptr;
        }

        inline double barrier_epsilon()
        {
            return general_problem_ptr->get_barrier_epsilon();
        }

        std::shared_ptr<IBarrierOptimizationSolver> inner_solver_ptr;
        std::unique_ptr<BarrierProblem> barrier_problem_ptr;
        OptimizationProblem* general_problem_ptr;
        int num_outer_iterations_;
        std::string name_;
    };

} // namespace opt
} // namespace ccd
