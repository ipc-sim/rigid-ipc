#include "ipc_solver.hpp"

#include <autodiff/autodiff_types.hpp>
#include <barrier/barrier.hpp>

namespace ccd {
namespace opt {

    IPCSolver::IPCSolver()
        : NewtonSolver()
    {
        convergence_criteria = ConvergenceCriteria::VELOCITY;
    }

    // Initialize the state of the solver using the settings saved in JSON
    void IPCSolver::settings(const nlohmann::json& json)
    {
        nlohmann::json edited_json = json;
        edited_json["convergence_criteria"] = "velocity";
        NewtonSolver::settings(edited_json);
        dhat_epsilon = json["dhat_epsilon"].get<double>();
    }

    // Export the state of the solver using the settings saved in JSON
    nlohmann::json IPCSolver::settings() const
    {
        nlohmann::json json = NewtonSolver::settings();
        json.erase("convergence_criteria");
        json["dhat_epsilon"] = dhat_epsilon;
        return json;
    }

    // Initialize the solver state for a new solve
    void IPCSolver::init_solve(const Eigen::VectorXd& x0)
    {
        assert(problem_ptr != nullptr);

        NewtonSolver::init_solve(x0);

        // Find a good initial value for κ
        // double min_barrier_stiffness = 1;
        double bbox_diagonal = problem_ptr->world_bbox_diagonal();
        double d0 = 1e-8 * bbox_diagonal;
        double min_barrier_stiffness = poly_log_barrier_hessian(d0, 1e-3);
        min_barrier_stiffness = 1.0e11 * /*average mass=*/1
            / (2e-8 * bbox_diagonal * min_barrier_stiffness);
        max_barrier_stiffness = 100 * min_barrier_stiffness;

        Eigen::VectorXd grad_B;
        int num_active_barriers;
        barrier_problem_ptr()->compute_barrier_term(
            x0, grad_B, num_active_barriers);
        double kappa;
        if (num_active_barriers > 0) {
            Eigen::VectorXd grad_E;
            barrier_problem_ptr()->compute_energy_term(x0, grad_E);

            kappa = -grad_B.dot(grad_E) / grad_B.squaredNorm();
            barrier_problem_ptr()->set_barrier_stiffness(std::min(
                max_barrier_stiffness, std::max(min_barrier_stiffness, kappa)));
        } else {
            barrier_problem_ptr()->set_barrier_stiffness(kappa = 1.0);
        }
        spdlog::info(
            "solver={} initial_num_active_barriers={:d} κ_min={:g} κ_max={:g} "
            "κ_g={:g} κ₀={:g}",
            name(), num_active_barriers, min_barrier_stiffness,
            max_barrier_stiffness, kappa,
            barrier_problem_ptr()->get_barrier_stiffness());
    }

    void IPCSolver::post_step_update()
    {
        // Adaptive κ
        double min_disti = problem_ptr->compute_min_distance(x_prev);
        double min_distj = problem_ptr->compute_min_distance(x);
        // Is the barrier having a difficulty pushing the bodies apart?
        if (min_disti < dhat_epsilon && min_distj < dhat_epsilon
            && min_distj < min_disti) {
            // Then increase the barrier stiffness.
            barrier_problem_ptr()->set_barrier_stiffness(std::min(
                max_barrier_stiffness,
                2 * barrier_problem_ptr()->get_barrier_stiffness()));
            spdlog::info(
                "solver={} iter={:d} msg=\"updated κ to {:g}\"", name(),
                iteration_number,
                barrier_problem_ptr()->get_barrier_stiffness());
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
