#include "split_distance_barrier_rb_problem.hpp"

#include <finitediff.hpp>

#include <constants.hpp>
#include <solvers/solver_factory.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

namespace opt {

    SplitDistanceBarrierRBProblem::SplitDistanceBarrierRBProblem()
        : DistanceBarrierRBProblem()
    {
    }

    // TODO: custom settings for split method

    ////////////////////////////////////////////////////////////
    // Rigid Body Problem

    void SplitDistanceBarrierRBProblem::simulation_step(
        bool& had_collision, bool& _has_intersections, bool solve_collision)
    {
        // Take an unconstrained time-step
        m_time_stepper->step(m_assembler, gravity, timestep());

        update_dof();

        had_collision =
            detect_collisions(poses_t0, poses_t1, CollisionCheck::CONSERVATIVE);

        // Check if minimum distance is violated
        min_distance = compute_min_distance(this->poses_to_dofs(poses_t1));
        if (min_distance >= 0) {
            spdlog::info("candidate_step min_distance={:.8e}", min_distance);

            // our constraint is really d > min_d, we want to run the
            // optimization when we end the step with small distances
            if (min_distance <= barrier_activation_distance()) {
                had_collision = true;
            }
        }

        if (solve_collision && had_collision) {
            update_constraint();
            opt::OptimizationResults result = solve_constraints();
            _has_intersections = take_step(result.x);
        } else {
            _has_intersections = has_intersections();
        }
    }

    ////////////////////////////////////////////////////////////
    // Barrier Problem

    // Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    double SplitDistanceBarrierRBProblem::compute_energy_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        Eigen::VectorXd diff = x - this->poses_to_dofs(poses_t1);
        Eigen::DiagonalMatrixXd M = m_assembler.m_rb_mass_matrix;

        if (compute_grad) {
            grad = M * diff;
#ifdef WITH_DERIVATIVE_CHECK
            Eigen::VectorXd grad_approx = eval_grad_energy_approx(*this, x);
            if (!fd::compare_gradient(grad, grad_approx)) {
                spdlog::error("finite gradient check failed for E(x)");
            }
#endif
        }

        if (compute_hess) {
            hess = Eigen::SparseDiagonal(M.diagonal());
#ifdef WITH_DERIVATIVE_CHECK
            Eigen::MatrixXd hess_approx = eval_hess_energy_approx(*this, x);
            if (!fd::compare_jacobian(hess, hess_approx)) {
                spdlog::error("finite hessian check failed for E(x)");
            }
#endif
        }

        return 0.5 * diff.transpose() * M * diff;
    }

} // namespace opt
} // namespace ccd
