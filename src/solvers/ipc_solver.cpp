#include "ipc_solver.hpp"

namespace ccd {
namespace opt {

    IPCSolver::IPCSolver()
        : NewtonSolver()
    {
    }

    // Initialize the state of the solver using the settings saved in JSON
    void IPCSolver::settings(const nlohmann::json& json)
    {
        NewtonSolver::settings(json);
        dhat_epsilon = json["dhat_epsilon"].get<double>();
    }

    // Export the state of the solver using the settings saved in JSON
    nlohmann::json IPCSolver::settings() const
    {
        nlohmann::json json = NewtonSolver::settings();
        json["dhat_epsilon"] = dhat_epsilon;
        return json;
    }

    // Initialize the solver state for a new solve
    void IPCSolver::init_solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);

        NewtonSolver::init_solve(x0);

        // Find a good initial value for κ
        double min_barrier_stiffness = 1; // TODO: Compute this value
        // Eigen::VectorX3d pointA = Eigen::VectorX3d::Zero(problem_ptr->dim());
        // Eigen::VectorX3d pointB = pointA;
        // pointB.y() += 1e-8 * bbox_diagonal;
        // min_barrier_stiffness =
        //     1e11 * barrier(point_point_distance(pointA,
        //     pointB)).getHessian();
        max_barrier_stiffness = 100 * min_barrier_stiffness;

        Eigen::VectorXd grad_B;
        int num_active_barriers;
        barrier_problem_ptr()->compute_barrier_term(
            x0, grad_B, num_active_barriers);
        if (num_active_barriers > 0) {
            Eigen::VectorXd grad_E;
            barrier_problem_ptr()->compute_energy_term(x0, grad_E);

            double kappa = -grad_B.dot(grad_E) / grad_B.squaredNorm();
            kappa = std::min(
                max_barrier_stiffness, std::max(min_barrier_stiffness, kappa));
            barrier_problem_ptr()->set_barrier_stiffness(kappa);
        } else {
            barrier_problem_ptr()->set_barrier_stiffness(1);
        }
    }

    void IPCSolver::post_step_update(
        const Eigen::VectorXd& xi, const Eigen::VectorXd& xj)
    {
        // Adaptive κ
        double min_disti = problem_ptr->compute_min_distance(xi);
        double min_distj = problem_ptr->compute_min_distance(xj);
        // Is the barrier having a difficulty pushing the bodies apart?
        if (min_disti < dhat_epsilon && min_distj < dhat_epsilon
            && min_distj < min_disti) {
            // Then increase the barrier stiffness.
            barrier_problem_ptr()->set_barrier_stiffness(std::min(
                max_barrier_stiffness,
                2 * barrier_problem_ptr()->get_barrier_stiffness()));
        }
    }

    // Solve the saved optimization problem to completion
    OptimizationResults IPCSolver::solve(const Eigen::VectorXd& x0)
    {
        init_solve(x0);
        return NewtonSolver::solve(x0);
    }

} // namespace opt
} // namespace ccd
