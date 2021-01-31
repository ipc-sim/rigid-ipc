// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include "homotopy_solver.hpp"

#include <barrier/barrier.hpp>
#include <utils/eigen_ext.hpp>

#include <solvers/solver_factory.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ipc::rigid {

HomotopySolver::HomotopySolver()
    : tinit(1.0)
    , t(1.0)
    , m(1.0)
    , c(0.01)
    , e_b(1e-5)
    , t_inc(2)
    , max_num_constraints(0)
    , num_outer_iterations(0)
{
    // inner solver to solve the unconstrained barrier problem
    inner_solver_ptr = std::make_unique<InnerNewtonSolver>();
}

void HomotopySolver::settings(const nlohmann::json& params)
{
    // termination criteria variables
    e_b = params["e_b"].get<double>();
    tinit = params["t_init"].get<double>();
    t_inc = params["t_inc"].get<double>();
    m = params["m"].get<double>();
    c = params["c"].get<double>();
}

nlohmann::json HomotopySolver::settings() const
{
    nlohmann::json json;
    json["e_b"] = e_b;
    json["t_inc"] = t_inc;
    json["t_init"] = tinit;
    json["m"] = m;
    json["c"] = c;
    // json["inner_solver"] = inner_solver_ptr->name();
    return json;
}

void HomotopySolver::init_solve(const Eigen::VectorXd& x0)
{
    assert(inner_solver_ptr != nullptr);
    assert(problem_ptr != nullptr);

    // Initialize the starting point as the starting point of the problem
    x0_i = x0;

    // Convert from the boolean vector to a vector of free dof indices
    inner_solver_ptr->init_solve(x0);

    // Reset the number of outer iterations
    num_outer_iterations = 0;

    // Find a good initial value for `t`
    Eigen::VectorXd grad_B;
    int num_active_barriers;
    problem_ptr->compute_barrier_term(x0_i, grad_B, num_active_barriers);
    if (false /*num_active_barriers > 0*/) {
        Eigen::VectorXd grad_E;
        problem_ptr->compute_energy_term(x0_i, grad_E);
        t = tinit = 1 / (-grad_B.dot(grad_E) / grad_B.squaredNorm());
    } else {
        t = tinit;
    }

    spdlog::debug(
        "solver={} d̂={} m={} t={} e_b={} c={} t_inc={}", name(),
        problem_ptr->barrier_activation_distance(), m, t, e_b, c, t_inc);
    problem_ptr->barrier_stiffness(1 / t);
}

OptimizationResults HomotopySolver::step_solve()
{
    assert(problem_ptr != nullptr);
    assert(inner_solver_ptr != nullptr);

    spdlog::debug(
        "solver={} it={} d̂={} m={} t={} e_b={} c={}", name(),
        num_outer_iterations, problem_ptr->barrier_activation_distance(), m, t,
        e_b, c);

    // propagate variables
    problem_ptr->barrier_stiffness(1 / t);
    inner_solver_ptr->e_b(e_b);
    inner_solver_ptr->c(c);
    inner_solver_ptr->t(t);
    inner_solver_ptr->m(m);

    OptimizationResults results = inner_solver_ptr->solve(x0_i);

    results.minf = problem_ptr->compute_energy_term(results.x);
    // Start next iteration from the ending optimal position
    x0_i = results.x;
    t *= t_inc;

    num_outer_iterations++;
    return results;
}

OptimizationResults HomotopySolver::solve(const Eigen::VectorXd& x0)
{
    init_solve(x0);

    OptimizationResults results;
    do {
        results = step_solve();
    } while (m / t > e_b);

    int num_constraints;
    problem_ptr->compute_barrier_term(x0_i, num_constraints);
    max_num_constraints = std::max(max_num_constraints, num_constraints);

    // make one last iteration with exactly eb
    t = m / e_b;
    const double t_used = t;
    results = step_solve();

    double min_dist = problem_ptr->compute_min_distance(results.x);
    spdlog::info(
        "solver={} t={:g} min_dist={:g} max_num_constraints={:d} {}", name(),
        t_used, min_dist, max_num_constraints,
        inner_solver_ptr->stats_string());

    spdlog::info(
        "solver={} c={:g} tinit={:g} tinc={:g} num_iterations={:d}", name(), c,
        tinit, t_inc, num_outer_iterations);

    return results;
}

} // namespace ipc::rigid
