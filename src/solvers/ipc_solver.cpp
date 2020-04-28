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
        num_kappa_updates = 0;
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
        double bbox_diagonal = problem_ptr->world_bbox_diagonal();
        double dhat = barrier_problem_ptr()->get_barrier_homotopy();
        double d0 = 1e-8 * bbox_diagonal;
        if (d0 >= dhat) {
            d0 = 0.5 * dhat; // TODO: this is untested
        }
        double min_barrier_stiffness =
            barrier_problem_ptr()->barrier_hessian(d0);
        min_barrier_stiffness = Constants::MIN_BARRIER_STIFFNESS_SCALE
            * problem_ptr->average_mass() / min_barrier_stiffness;
        if (!std::isfinite(min_barrier_stiffness)) {
            spdlog::error(
                "κ_min is not finite using a distance d₀={:g} with d̂={:g}", d0,
                dhat);
            throw NotImplementedError("κ_min is not finite");
        }
#ifdef USE_DISTANCE_SQUARED
        min_barrier_stiffness /= 4 * d0 * d0;
#endif
        max_barrier_stiffness = 100 * min_barrier_stiffness;

        Eigen::VectorXd grad_B;
        int num_active_barriers;
        barrier_problem_ptr()->compute_barrier_term(
            x0, grad_B, num_active_barriers);
        double kappa = 1.0;
        if (num_active_barriers > 0 && grad_B.squaredNorm() > 0) {
            Eigen::VectorXd grad_E;
            barrier_problem_ptr()->compute_energy_term(x0, grad_E);
            kappa = -grad_B.dot(grad_E) / grad_B.squaredNorm();
            assert(std::isfinite(kappa));
        }
        barrier_problem_ptr()->set_barrier_stiffness(std::min(
            max_barrier_stiffness, std::max(min_barrier_stiffness, kappa)));
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
            num_kappa_updates++;
        }
    }

    // Solve the saved optimization problem to completion
    OptimizationResults IPCSolver::solve(const Eigen::VectorXd& x0)
    {
        init_solve(x0);
        OptimizationResults results = NewtonSolver::solve(x0);
        std::cout
            << fmt::format(
                   "GREP_ME solver={} min_dist={:g} num_kappa_updates={:d} {}",
                   name(), problem_ptr->compute_min_distance(results.x),
                   num_kappa_updates, NewtonSolver::stats())
            << std::endl;
        return results;
    }

} // namespace opt
} // namespace ccd
