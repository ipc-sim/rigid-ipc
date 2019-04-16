#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/lcp_solver.hpp>

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
        const Eigen::VectorXd& b,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>& g,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd& x)>& jac_g,
        const int max_iter, const callback_intermediate_ncp& callback,
        const NcpUpdate update_type, const LCPSolver lcp_solver,
        Eigen::VectorXd& x, Eigen::VectorXd& alpha,
        const bool check_convergence, const bool check_convergence_unfeasible,
        const double convergence_tolerance = 1E-10);

} // namespace opt
} // namespace ccd
