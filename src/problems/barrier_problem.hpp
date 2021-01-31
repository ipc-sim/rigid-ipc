#pragma once

#include <memory>

#include <solvers/optimization_solver.hpp>

namespace ipc::rigid {

class BarrierProblem : public virtual OptimizationProblem {
public:
    virtual ~BarrierProblem() = default;

    /// Compute the objective function f(x)
    virtual double compute_objective(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true) override;

    // Include thes lines to avoid issues with overriding inherited
    // functions with the same name.
    // (http://www.cplusplus.com/forum/beginner/24978/)
    using OptimizationProblem::compute_objective;

    /// @brief Compute \f$E(x)\f$ in
    /// \f$f(x) = E(x) + \kappa \sum_{k \in C} b(d(x_k))\f$
    /// @param[in]  x     input value
    /// @param[out] grad  \f$\nabla E(x)\f$ to compute
    /// @param[out] hess  \f$\nabla^2 E(x)\f$ to compute
    /// @param[in]  compute_grad  True if \f$\nabla E(x)\f$ should be
    ///                           computed
    /// @param[in]  compute_hess  True if \f$\nabla^2 E(x)\f$ should be
    ///                           computed
    virtual double compute_energy_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true) = 0;

    /// @brief Compute \f$\sum_{k ∈ C} b(d(x_k))\f$ in
    ///        \f$f(x) = E(x) + κ \sum_{k \in C} b(d(x_k))\f$
    /// @param[in]  x     input value
    /// @param[out] grad  \f$\nabla \sum_{k ∈ C} b(d(x_k))\f$ to compute
    /// @param[out] hess  \f$\nabla^2 \sum_{k ∈ C} b(d(x_k))\f$ to compute
    /// @param[in]  compute_grad  True if \f$\nabla \sum_{k ∈ C}
    ///                           b(d(x_k))\f$ should be computed
    /// @param[in]  compute_hess  True if \f$\nabla^2 \sum_{k ∈ C}
    ///                           b(d(x_k))\f$ should be computed
    /// @returns number of active barriers
    virtual double compute_barrier_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        int& num_constraints,
        bool compute_grad = true,
        bool compute_hess = true) = 0;

    // --------------------------------------------------------------------
    // Convience functions
    // --------------------------------------------------------------------

    /// @brief Compute \f$E(x)\f$ in
    /// \f$f(x) = E(x) + \kappa \sum_{k \in C} b(d(x_k))\f$
    ///
    /// Does not compute the gradient or hessian.
    ///
    /// @param[in]  x     input value
    virtual double compute_energy_term(const Eigen::VectorXd& x) final
    {
        Eigen::VectorXd grad;
        Eigen::SparseMatrix<double> hess;
        return compute_energy_term(
            x, grad, hess,
            /*compute_grad=*/false, /*compute_hess=*/false);
    }

    /// @brief Compute \f$E(x)\f$ in
    /// \f$f(x) = E(x) + \kappa \sum_{k \in C} b(d(x_k))\f$
    ///
    /// Does not compute the hessian.
    ///
    /// @param[in]  x     input value
    /// @param[out] grad  \f$\nabla E(x)\f$ to compute
    virtual double
    compute_energy_term(const Eigen::VectorXd& x, Eigen::VectorXd& grad) final
    {
        Eigen::SparseMatrix<double> hess;
        return compute_energy_term(
            x, grad, hess,
            /*compute_grad=*/true, /*compute_hess=*/false);
    }

    /// @brief Compute \f$\sum_{k ∈ C} b(d(x_k))\f$ in
    ///        \f$f(x) = E(x) + κ \sum_{k \in C} b(d(x_k))\f$
    ///
    /// Does not compute the gradient or hessian.
    ///
    /// @param[in]  x     input value
    /// @returns number of active barriers
    virtual double
    compute_barrier_term(const Eigen::VectorXd& x, int& num_constraints) final
    {
        Eigen::VectorXd grad;
        Eigen::SparseMatrix<double> hess;
        return compute_barrier_term(
            x, grad, hess, num_constraints,
            /*compute_grad=*/false, /*compute_hess=*/false);
    }

    /// @brief Compute \f$\sum_{k ∈ C} b(d(x_k))\f$ in
    ///        \f$f(x) = E(x) + κ \sum_{k \in C} b(d(x_k))\f$
    ///
    /// Does not compute the hessian.
    ///
    /// @param[in]  x     input value
    /// @param[out] grad  \f$\nabla \sum_{k ∈ C} b(d(x_k))\f$ to compute
    /// @returns number of active barriers
    virtual double compute_barrier_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        int& num_constraints) final
    {
        Eigen::SparseMatrix<double> hess;
        return compute_barrier_term(
            x, grad, hess, num_constraints,
            /*compute_grad=*/true, /*compute_hess=*/false);
    }

    virtual bool is_barrier_problem() const override { return true; }

    /// Compute the value of the barrier hessian at a value x
    virtual double barrier_hessian(double x) const = 0;

    virtual double barrier_activation_distance() const = 0;
    virtual void barrier_activation_distance(const double dhat) = 0;

    virtual double barrier_stiffness() const = 0;
    virtual void barrier_stiffness(const double kappa) = 0;

protected:
    bool m_use_barriers = true;
};

/// Helper Functions for checking finite differences
/// ------------------------------------------------
Eigen::VectorXd
eval_grad_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x);
Eigen::MatrixXd
eval_hess_energy_approx(BarrierProblem& problem, const Eigen::VectorXd& x);

// WARNING: You cannot create a eval_*_barrier_approx because barrier term
// eval re-computes the constraint set which can lead to poor performance
// and discontiuties when using finite differences.

} // namespace ipc::rigid
