#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/OptimizationProblem.hpp>
#include <opt/OptimizationResults.hpp>
#include <opt/lcp_solver.hpp>

namespace ccd {
namespace opt {

    typedef std::function<void(const Eigen::VectorXd& x,
        const Eigen::VectorXd& alpha, const double gamma)>
        callback_intermediate_ncp;

    enum UpdateType {
        GRADIENT_ONLY,
        OTHER,
    };

    /**
     * @brief solves the optimization problem
     *      Ax = b + \nabla g(x)^T \alpha
     *      alpha >=0, g(x) >=0   alpha^T g(x) = 0
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

     */
    bool solve_ncp(const Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
        const callback_g& g, const callback_jac_g& jac_g, const int max_iter,
        const callback_intermediate_ncp& callback, const UpdateType update_type,
        const LCPSolver lcp_solver, Eigen::VectorXd& x, Eigen::VectorXd& alpha,
        const bool check_convergence);

} // namespace opt
} // namespace ccd
