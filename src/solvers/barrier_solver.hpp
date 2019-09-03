#pragma once

/**
 * @brief Solve the optimization problem using Newton's Method with barriers for
 * the constraints.
 */
#include <memory>

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

    class BarrierProblem : public virtual IBarrierProblem {
    public:
        BarrierProblem(IBarrierGeneralProblem& problem);
        ~BarrierProblem() override = default;

        double eval_f(const Eigen::VectorXd& x) override;

        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;

        double get_barrier_epsilon() override
        {
            return general_problem->get_barrier_epsilon();
        }

        const Eigen::VectorXd& starting_point() override { return x0; }
        const int& num_vars() override { return num_vars_; }
        const Eigen::VectorXb& is_dof_fixed() override;

        Multiprecision eval_mp_f(const Eigen::VectorXd& x) override;

        bool has_collisions(const Eigen::VectorXd& sigma_i,
            const Eigen::VectorXd& sigma_j) const override;

        // new functions for termiantion criteria
        double get_termination_threshold() const override;
        Eigen::VectorXd eval_grad_E(const Eigen::VectorXd& xk) override;
        Eigen::VectorXd eval_grad_B(const Eigen::VectorXd& xk, int&) override;

#ifdef DEBUG_LINESEARCH
        Eigen::MatrixXi debug_edges() const override
        {
            return general_problem->debug_edges();
        }
        Eigen::MatrixXd debug_vertices(
            const Eigen::VectorXd& sigma) const override
        {
            return general_problem->debug_vertices(sigma);
        }
        Eigen::MatrixXd debug_vertices_t0() const override
        {
            return general_problem->debug_vertices_t0();
        }

        double debug_min_distance(const Eigen::VectorXd& sigma) const override
        {
            return general_problem->debug_min_distance(sigma);
        }
#endif

        IBarrierGeneralProblem* general_problem;
        Eigen::VectorXd x0;

        int num_vars_;
        double inner_solver_threshold;

        double t;

    };

    class BarrierSolver : public virtual IStateOptimizationSolver {
    public:
        BarrierSolver();
        BarrierSolver(const std::string& name);
        ~BarrierSolver() override {}

        void set_problem(IBarrierGeneralProblem& problem);

        OptimizationResults solve() override;
        void init_solve() override;
        OptimizationResults step_solve() override;

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        const std::string& name() const override { return name_; }
        bool has_inner_solver() override { return true; }
        const IBarrierOptimizationSolver& inner_solver() override
        {
            return *inner_solver_ptr;
        }
        void inner_solver_settings(const nlohmann::json& json);

        Eigen::VectorXd get_grad_kkt() const override;
        int num_outer_iterations() const override
        {
            return num_outer_iterations_;
        }

        double min_barrier_epsilon;
        int max_iterations;
        double tinit;
        double t;
        double m;
        double e_b;
        double t_inc;

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
        IBarrierGeneralProblem* general_problem_ptr;
        int num_outer_iterations_;
        std::string name_;

        std::stringstream debug;
    };

} // namespace opt
} // namespace ccd
