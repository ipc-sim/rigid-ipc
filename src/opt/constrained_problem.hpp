#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

#include <opt/optimization_problem.hpp>

namespace ipc::rigid {

class ConstrainedProblem : public virtual OptimizationProblem {
public:
    virtual ~ConstrainedProblem() = default;

    /// Compute the constraint function g(x) ≥ 0
    virtual void compute_constraints(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& jac_gx,
        std::vector<Eigen::SparseMatrix<double>>& hess_gx,
        bool compute_grad = true,
        bool compute_hess = true) = 0;

    /// Compute the constraint function g(x) ≥ 0
    virtual void compute_constraints_using_normals(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& jac_gx) = 0;

    // --------------------------------------------------------------------
    // Convience functions
    // --------------------------------------------------------------------

    virtual void
    compute_constraints(const Eigen::VectorXd& x, Eigen::VectorXd& gx) final
    {
        Eigen::MatrixXd jac_gx;
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        return compute_constraints(
            x, gx, jac_gx, hess_gx,
            /*compute_grad=*/false, /*compute_hess=*/false);
    }

    virtual void compute_constraints(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& jac_gx) final
    {
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        return compute_constraints(
            x, gx, jac_gx, hess_gx,
            /*compute_grad=*/true, /*compute_hess=*/false);
    }

    virtual bool is_constrained_problem() const { return true; }
};

} // namespace ipc::rigid
