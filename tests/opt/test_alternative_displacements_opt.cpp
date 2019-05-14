#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <state.hpp>

#include <Eigen/Geometry>

#define NUM_ANGLES (5)

using namespace ccd;
using namespace opt;

TEST_CASE("test the setup", "[opt][displacements][barrier]")
{
    State state;
    state.constraint_function
        = GENERATE(ccd::ConstraintType::BARRIER, ccd::ConstraintType::VOLUME);

    state.vertices.resize(4, 2);
    state.displacements.resize(4, 2);
    state.edges.resize(2, 2);

    state.vertices << -1, 1, 1, 1, -2, 0, 2, 0;
    state.displacements << 0, -2, 0, -2, 0, 0, 0, 0;
    state.edges << 0, 1, 2, 3;

    state.recompute_collision_set = GENERATE(false, true);
    state.detection_method = DetectionMethod::BRUTE_FORCE;
    state.solver_settings.barrier_epsilon = 1.0;
    state.solver_settings.method = OptimizationMethod::BARRIER_NEWTON;

    state.reset_optimization_problem();

    CHECK(state.opt_problem.num_vars == state.displacements.size());
    CHECK(state.opt_problem.num_constraints == state.edges.rows());

    CHECK(state.opt_problem.x0.size() == state.opt_problem.num_vars);
    CHECK(state.opt_problem.x_lower.size() == state.opt_problem.num_vars);
    CHECK(state.opt_problem.x_upper.size() == state.opt_problem.num_vars);
    CHECK(
        state.opt_problem.g_lower.size() == state.opt_problem.num_constraints);
    CHECK(
        state.opt_problem.g_upper.size() == state.opt_problem.num_constraints);

    state.opt_problem.x0.setRandom();

    // Test ∇f
    Eigen::VectorXd finite_grad(state.opt_problem.num_vars);
    finite_gradient(state.opt_problem.x0, state.opt_problem.func_f(), finite_grad);
    Eigen::VectorXd analytic_grad
        = state.opt_problem.eval_grad_f(state.opt_problem.x0);
    CHECK(compare_gradient(finite_grad, analytic_grad));

    // Test ∇²f
    Eigen::MatrixXd finite_hessian(
        state.opt_problem.num_vars, state.opt_problem.num_vars);
    finite_jacobian(
        state.opt_problem.x0, state.opt_problem.func_grad_f(), finite_hessian);
    Eigen::MatrixXd analytic_hessian
        = state.opt_problem.eval_hessian_f(state.opt_problem.x0);
    CHECK(compare_jacobian(finite_hessian, analytic_hessian));

    // Test ∇g
    Eigen::MatrixXd finite_jac(
        state.opt_problem.num_constraints, state.opt_problem.num_vars);
    finite_jacobian(state.opt_problem.x0, state.opt_problem.func_g(), finite_jac);
    Eigen::MatrixXd analytic_jac
        = state.opt_problem.eval_jac_g(state.opt_problem.x0);
    CHECK(compare_jacobian(finite_jac, analytic_jac));

    // Test ∇²g
    Eigen::MatrixXd finite_hessian_gi(
        state.opt_problem.num_vars, state.opt_problem.num_vars);
    finite_hessian_gi.setOnes();
    long i;
    auto diff_i = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        auto jac = state.opt_problem.eval_jac_g(x);
        return jac.row(i);
    };
    std::vector<Eigen::MatrixXd> analytic_hessian_g
        = state.opt_problem.eval_hessian_g(state.opt_problem.x0);
    for (i = 0; i < long(analytic_hessian_g.size()); i++) {
        finite_jacobian(state.opt_problem.x0, diff_i, finite_hessian_gi);
        CHECK(compare_jacobian(
            finite_hessian_gi, analytic_hessian_g[unsigned(i)]));
    }
}

TEST_CASE("two rotating edges", "[opt][displacements][barrier]")
{
    State state;
    state.load_scene(std::string(FIXTURES_DIR) + "/single-falling-edge.json");

    REQUIRE(state.vertices.rows() == 4);
    REQUIRE(state.displacements.rows() == 4);
    REQUIRE(state.edges.rows() == 2);

    state.constraint_function = ccd::ConstraintType::BARRIER;

    state.recompute_collision_set = GENERATE(false, true);
    state.detection_method = DetectionMethod::BRUTE_FORCE;
    state.solver_settings.method = OptimizationMethod::BARRIER_NEWTON;

    double theta1 = 2 * M_PI / NUM_ANGLES * GENERATE(range(0, NUM_ANGLES));
    double theta2 = 2 * M_PI / NUM_ANGLES * GENERATE(range(0, NUM_ANGLES));

    for (int i = 0; i < state.edges.rows(); i++) {
        Eigen::Rotation2D<double> R(i == 0 ? theta1 : theta2);
        Eigen::RowVector2d center = (state.vertices.row(state.edges(i, 1))
                                        + state.vertices.row(state.edges(i, 0)))
            / 2;
        state.vertices.row(state.edges(i, 0))
            = (R * (state.vertices.row(state.edges(i, 0)) - center).transpose())
                  .transpose()
            + center;
        state.vertices.row(state.edges(i, 1))
            = (R * (state.vertices.row(state.edges(i, 1)) - center).transpose())
                  .transpose()
            + center;
    }

    state.solver_settings.min_barrier_epsilon = 1e-3;
    state.solver_settings.line_search_tolerance = 1e-4;
    state.solver_settings.absolute_tolerance = 1e-3;

    state.optimize_displacements("");
    CHECK(state.opt_results.success);
}

TEST_CASE("corner case", "[opt][displacements][barrier]")
{
    State state;
    state.load_scene(std::string(FIXTURES_DIR) + "/corner-case.json");

    REQUIRE(state.vertices.rows() == 5);
    REQUIRE(state.displacements.rows() == 5);
    REQUIRE(state.edges.rows() == 3);

    state.constraint_function = ccd::ConstraintType::BARRIER;

    state.recompute_collision_set = GENERATE(false, true);
    state.detection_method = DetectionMethod::BRUTE_FORCE;
    state.solver_settings.method = OptimizationMethod::BARRIER_NEWTON;

    state.solver_settings.min_barrier_epsilon = 1e-3;
    state.solver_settings.line_search_tolerance = 1e-4;
    state.solver_settings.absolute_tolerance = 1e-3;

    state.optimize_displacements("");
    CHECK(state.opt_results.success);
}
