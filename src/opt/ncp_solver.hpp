#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/lcp_solver.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>
#include <opt/solver.hpp>

namespace ccd {
namespace opt {

    typedef std::function<void(const Eigen::VectorXd& x,
        const Eigen::VectorXd& alpha, const double gamma)>
        callback_intermediate_ncp;

    /**
     * @brief Method for updating candidate solution during NCP iterations
     *
     * Update: \f$ x_{i+1} = x_{i} + \Delta x \f$
     */
    enum class NcpUpdate {
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i\f$
        G_GRADIENT,
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i + A^{-1}b - x_i\f$
        LINEARIZED
    };

    static const char* NcpUpdateNames[] = { "G_GRADIENT", "LINEARIZED" };

    class NCPSolver : public OptimizationSolver{
    public:
        NCPSolver();
        ~NCPSolver() override;
        NCPSolver(const bool keep_in_unfeasible, const bool check_convergence,
            const double convergence_tolerance, const NcpUpdate update_type,
            const LCPSolver lcp_solver, const int max_iterations);

        OptimizationResults solve(OptimizationProblem& problem) override;
        bool solve_ncp(const Eigen::SparseMatrix<double>& A,
            const Eigen::VectorXd& b, OptimizationProblem& problem,
            Eigen::VectorXd& x_opt, Eigen::VectorXd& alpha_opt);
        void compute_linear_system(OptimizationProblem& problem);

        bool compute();
        void initialize();
        void solve_lcp(Eigen::VectorXd& delta_x);
        void move_to_unfeasible_domain(Eigen::VectorXd& delta_x, double& gamma);
        void update_candidate(Eigen::VectorXd& delta_x, double& gamma);

        Eigen::VectorXd Ainv(const Eigen::VectorXd& x);

        // ---------------------
        // Configuration
        // ---------------------
        bool keep_in_unfeasible;
        bool check_convergence;
        double convergence_tolerance;
        NcpUpdate update_type;
        LCPSolver lcp_solver;
        int max_iterations;

        // ---------------------
        // Optimization Specifics
        // ---------------------
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        OptimizationProblem* problem;
        //        callback_intermediate_ncp& callback;

        // -----------------------
        // Optimization Status
        // -----------------------
        std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> Asolver;
        Eigen::VectorXd g_xi;
        Eigen::MatrixXd jac_g_xi;

        // ----------------------
        // Optimization results
        // ----------------------
        Eigen::VectorXd xi;
        Eigen::VectorXd alpha_i;
    };

    /**
     * @brief solves the optimization problem
     *      Ax = b + \nabla g(x)^T \alpha           (1)
     *      alpha >=0, g(x) >=0   alpha^T g(x) = 0  (2)
     *
     * It can be used to solve the optimization problem of the form
     *      Min 1/2 || Ax + b ||^2
     *      s.t     g(x) >=0
     *
     * @param A         : Sparse Square Matrix (dof x dof)
     * @param b         : Vector (dof)
     * @param g         : constraint function g(x)
     * @param jac_g     : jacobian of constraint function \nabla g(x)
     * @param max_iter  : max-iterations for the newton-like loop
     * @param callback  : callback function called after each iteration
     * @param x         : the solution of the problem.
     * @param alpha     : the solution of the dual problem.
     *
     * @param check_convergence  :
     *                  Solver ends when (1) and (2) are satisfied
     * @param check_convergence_unfeasible:
     *                  Used for g(x) that well defined only for g(x) <=0.
     *                  Solver ends when (1) is satisfied but we stay in the
     *                  unfeasible domain until an x satisfied (2)
     * @param convergence_tolerance: used to evaluate (2)
     */
    bool solve_ncp(const Eigen::SparseMatrix<double>& A,
        const Eigen::VectorXd& b, OptimizationProblem& problem,
        const int max_iter, const callback_intermediate_ncp& callback,
        const NcpUpdate update_type, const LCPSolver lcp_solver,
        Eigen::VectorXd& x, Eigen::VectorXd& alpha,
        const bool keep_in_unfeasible, const bool check_convergence,
        const double convergence_tolerance = 1E-10);

} // namespace opt
} // namespace ccd
