#include "ipc_solver.hpp"

#include <ipc/barrier/adaptive_stiffness.hpp>

#include <autodiff/autodiff_types.hpp>
#include <barrier/barrier.hpp>

namespace ipc::rigid {

IPCSolver::IPCSolver()
    : NewtonSolver()
{
    convergence_criteria = ConvergenceCriteria::VELOCITY;
}

// Initialize the state of the solver using the settings saved in JSON
void IPCSolver::settings(const nlohmann::json& json)
{
    NewtonSolver::settings(json);
    dhat_epsilon = json["dhat_epsilon"].get<double>();
    min_barrier_stiffness_scale =
        json["min_barrier_stiffness_scale"].get<double>();
    num_kappa_updates = 0;
}

// Export the state of the solver using the settings saved in JSON
nlohmann::json IPCSolver::settings() const
{
    nlohmann::json json = NewtonSolver::settings();
    json["dhat_epsilon"] = dhat_epsilon;
    json["min_barrier_stiffness_scale"] = min_barrier_stiffness_scale;
    return json;
}

// Initialize the solver state for a new solve
void IPCSolver::init_solve(const Eigen::VectorXd& x0)
{
    assert(problem_ptr != nullptr);

    NewtonSolver::init_solve(x0);

    double bbox_diagonal = problem_ptr->world_bbox_diagonal();
    double dhat = barrier_problem_ptr()->barrier_activation_distance();
    double average_mass = problem_ptr->average_mass();

    Eigen::VectorXd grad_E;
    barrier_problem_ptr()->compute_energy_term(x0, grad_E);

    Eigen::VectorXd grad_B;
    int num_active_barriers;
    barrier_problem_ptr()->compute_barrier_term(
        x0, grad_B, num_active_barriers);

    double kappa = initial_barrier_stiffness(
        bbox_diagonal, dhat, average_mass, grad_E, grad_B,
        max_barrier_stiffness, min_barrier_stiffness_scale);

    barrier_problem_ptr()->barrier_stiffness(kappa);
    spdlog::info(
        "solver={} initial_num_active_barriers={:d} κ₀={:g} κ_max={:g}", name(),
        num_active_barriers, kappa, max_barrier_stiffness);

    // Compute the initial minimum distance
    prev_min_distance = problem_ptr->compute_min_distance(x0);
}

void IPCSolver::post_step_update()
{
    NewtonSolver::post_step_update();

    // Adaptive κ
    double min_distance = problem_ptr->compute_min_distance(x);
    spdlog::debug(
        "solver={} iter={:d} min_distance={:g}", name(), iteration_number,
        min_distance);
    double kappa = barrier_problem_ptr()->barrier_stiffness();
    update_barrier_stiffness(
        prev_min_distance, min_distance, max_barrier_stiffness, kappa,
        problem_ptr->world_bbox_diagonal(), dhat_epsilon);
    if (kappa != barrier_problem_ptr()->barrier_stiffness()) {
        spdlog::info(
            "solver={} iter={:d} msg=\"updated κ to {:g}\"", name(),
            iteration_number, kappa);
        barrier_problem_ptr()->barrier_stiffness(kappa);
        num_kappa_updates++;
    }
    prev_min_distance = min_distance;
}

// Solve the saved optimization problem to completion
OptimizationResults IPCSolver::solve(const Eigen::VectorXd& x0)
{
    init_solve(x0);
    OptimizationResults results = NewtonSolver::solve(x0);
    spdlog::info(
        "solver={} min_dist={:g} {}", name(),
        problem_ptr->compute_min_distance(results.x), stats_string());
    return results;
}

std::string IPCSolver::stats_string() const
{
    return fmt::format(
        "num_kappa_updates={:d} {}", num_kappa_updates,
        NewtonSolver::stats_string());
}

nlohmann::json IPCSolver::stats() const
{
    nlohmann::json stats_json = NewtonSolver::stats();
    stats_json["num_kappa_updates"] = num_kappa_updates;
    return stats_json;
}

} // namespace ipc::rigid
