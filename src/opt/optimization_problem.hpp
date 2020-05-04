#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    class OptimizationProblem {
    public:
        virtual ~OptimizationProblem() = default;

        /// Compute the objective function f(x)
        virtual double compute_objective(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad = true,
            bool compute_hess = true) = 0;

        // --------------------------------------------------------------------
        // Convience functions
        // --------------------------------------------------------------------

        virtual double compute_objective(const Eigen::VectorXd& x) final
        {
            Eigen::VectorXd grad;
            Eigen::SparseMatrix<double> hess;
            return compute_objective(
                x, grad, hess,
                /*compute_grad=*/false, /*compute_hess=*/false);
        }

        virtual double
        compute_objective(const Eigen::VectorXd& x, Eigen::VectorXd& grad) final
        {
            Eigen::SparseMatrix<double> hess;
            return compute_objective(
                x, grad, hess,
                /*compute_grad=*/true, /*compute_hess=*/false);
        }

        /// @returns the number of variables
        virtual int num_vars() const = 0;

        /// @returns A vector of booleans indicating if a DoF is fixed.
        virtual const Eigen::VectorXb& is_dof_fixed() = 0;

        /// Determine if there is a collision between two configurations
        virtual bool has_collisions(
            const Eigen::VectorXd& xi, const Eigen::VectorXd& xj) const = 0;

        /// Compute the minimum distance among geometry
        virtual double compute_min_distance(const Eigen::VectorXd& x) const = 0;

        /// Get the world coordinates of the vertices
        virtual Eigen::MatrixXd
        world_vertices(const Eigen::VectorXd& x) const = 0;

        /// Get the length of the diagonal of the worlds bounding box
        virtual double world_bbox_diagonal() const = 0;

        /// Get the average mass
        virtual double average_mass() const = 0;

        /// Get the time-step
        virtual double timestep() const = 0;

        virtual bool is_barrier_problem() const { return false; }
        virtual bool is_constrained_problem() const { return false; }
    };

    /// Helper Functions for checking finite differences
    /// ------------------------------------------------
    Eigen::VectorXd eval_grad_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x);
    Eigen::MatrixXd eval_hess_objective_approx(
        OptimizationProblem& problem, const Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
