#pragma once
#include <memory> // shared_ptr

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/constrained_problem.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <solvers/lcp_solver.hpp>
#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {

    /**
     * @brief Method for updating candidate solution during NCP iterations
     *
     * Update: \f$ x_{i+1} = x_{i} + \Delta x \f$
     */
    enum class NCPUpdate {
        /// \f$\Delta x = A^{-1} \nabla g(x_i)^T \lambda_i\f$
        G_GRADIENT,
        /// \f$\Delta x = A^{-1} \nabla g(x_i)^T \lambda_i + A^{-1}b - x_i\f$
        LINEAR
    };

    class NCPSolver : public OptimizationSolver {
    public:
        NCPSolver();
        virtual ~NCPSolver() = default;

        /// Initialize the state of the solver using the settings saved in JSON
        virtual void settings(const nlohmann::json& json) override;
        /// Export the state of the solver using the settings saved in JSON
        virtual nlohmann::json settings() const override;

        /// An identifier for the solver class
        static std::string solver_name() { return "ncp_solver"; }
        /// An identifier for this solver
        virtual std::string name() const override
        {
            return NCPSolver::solver_name();
        }

        void set_problem(OptimizationProblem& problem) override
        {
            assert(problem.is_constrained_problem());
            problem_ptr = dynamic_cast<ConstrainedProblem*>(&problem);
        }

        /// Initialize the solver state for a new solve
        void init_solve(const Eigen::VectorXd& x0) override;
        /// Solve the saved optimization problem to completion
        virtual OptimizationResults solve(const Eigen::VectorXd& x0) override;
        OptimizationResults
        solve(const Eigen::VectorXd& x0, const bool use_grad);
        /// Perform a single step of solving the optimization problem
        virtual OptimizationResults step_solve() override;

        // --------------------------------------------------------------------
        // Configuration
        // --------------------------------------------------------------------
        bool do_line_search;
        bool solve_for_active_cstr;
        double convergence_tolerance;
        NCPUpdate update_type;
        LCPSolver lcp_solver;
        int max_iterations;

        // --------------------------------------------------------------------
        // Optimization Status
        // --------------------------------------------------------------------
        std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> Asolver;
        Eigen::VectorXd g_xi;
        Eigen::MatrixXd jac_g_xi;

        // --------------------------------------------------------------------
        // Optimization results
        // --------------------------------------------------------------------
        Eigen::VectorXd xi;
        Eigen::VectorXd lambda_i;

    protected:
        void compute_linear_system(ConstrainedProblem& problem_ptr_);
        void compute_initial_solution();

        /**
         * @brief Linearize the problem and solve for primal variables (xᵢ₊₁)
         * and dual variables (λᵢ).
         *
         * Linearization:
         * \f{aligned}{
         *      A x_{i+1} = b + \nabla g(x_i)^T \lambda_i \\
         *      0 \leq \lambda_i \perp g(x_i) + \nabla g(x_i) Δx \geq 0
         * \f}
         * Update:
         * \f{aligned}{
         *      x_{i+1} = x_i + \Delta x
         * \f}
         * \f$\Delta x\f$:
         * * g_gradient update: \f$\Delta x = A^{-1} [\nabla g(x_i)]^T \lambda_i
         *      \f$
         * * linearized update: \f$\Delta x = A^{-1} [\nabla g(x_i)]^T \lambda_i
         *      + A^{-1}b - x_i\f$
         *
         * We want to take our problem to the form
         * \f{aligned}{
         *      s = q + N (M\lambda_i + p) \\
         *      0 \leq \lambda_i \perp s \geq 0
         * \f}
         * where
         * \f{aligned}{
         *      q &= g(x_i) \\
         *      N &= \nabla g(x_i) \\
         *      M &= A^{-1} [\nabla g(x_i)]^{T} \\
         *      p &= \Delta x - A^{-1} [\nabla g(x_i)]^{T} \lambda_i
         * \f}
         */
        Eigen::VectorXd solve_lcp();

        /**
         * @brief Compute \f$A^{-1}x\f$.
         *
         * Uses a precomputed decomposition of A to solve the system in
         * \f$O(n^2)\f$.
         */
        Eigen::VectorXd Ainv(const Eigen::VectorXd& x) const;

        // --------------------------------------------------------------------
        // Fields
        // --------------------------------------------------------------------
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        ConstrainedProblem* problem_ptr;

        /// @brief Current number of outer iterations completed.
        int num_outer_iterations_;
        std::string name_;
        bool m_use_gradient = true;
    };

    void
    zero_out_fixed_dof(const Eigen::VectorXb& is_fixed, Eigen::MatrixXd& jac);
} // namespace opt
} // namespace ccd
