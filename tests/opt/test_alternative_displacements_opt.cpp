#include <iostream>

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <state.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Test the setup", "[opt][displacements]")
{
    State state;
    state.use_alternative_formulation = GENERATE(false, true);

    state.vertices.resize(4, 2);
    state.displacements.resize(4, 2);
    state.edges.resize(2, 2);

    state.vertices << -1, 1, 1, 1, -2, 0, 2, 0;
    state.displacements << 0, -2, 0, -2, 0, 0, 0, 0;
    state.edges << 0, 1, 2, 3;

    state.recompute_collision_set = GENERATE(false, true);
    state.detection_method = DetectionMethod::BRUTE_FORCE;

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
    finite_gradient(state.opt_problem.x0, state.opt_problem.f, finite_grad);
    Eigen::VectorXd analytic_grad
        = state.opt_problem.grad_f(state.opt_problem.x0);
    CHECK(compare_gradient(finite_grad, analytic_grad));

    // Test ∇²f
    Eigen::MatrixXd finite_hessian(
        state.opt_problem.num_vars, state.opt_problem.num_vars);
    finite_jacobian(
        state.opt_problem.x0, state.opt_problem.grad_f, finite_hessian);
    Eigen::MatrixXd analytic_hessian
        = state.opt_problem.hessian_f(state.opt_problem.x0);
    CHECK(compare_jacobian(finite_hessian, analytic_hessian));

    // Test ∇g
    Eigen::MatrixXd finite_jac(
        state.opt_problem.num_constraints, state.opt_problem.num_vars);
    finite_jacobian(state.opt_problem.x0, state.opt_problem.g, finite_jac);
    Eigen::MatrixXd analytic_jac
        = state.opt_problem.jac_g(state.opt_problem.x0);
    CHECK(compare_jacobian(finite_jac, analytic_jac));

    // Test ∇²g
    Eigen::MatrixXd finite_hessian_gi(
        state.opt_problem.num_vars, state.opt_problem.num_vars);
    finite_hessian_gi.setOnes();
    long i;
    auto diff_i = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        assert(i >= 0 && i < state.opt_problem.num_constraints);
        auto jac = state.opt_problem.jac_g(x);
        return jac.row(i);
    };
    std::vector<Eigen::MatrixXd> analytic_hessian_g
        = state.opt_problem.hessian_g(state.opt_problem.x0);
    for (i = 0; i < state.opt_problem.num_constraints; i++) {
        finite_jacobian(state.opt_problem.x0, diff_i, finite_hessian_gi);
        CHECK(compare_jacobian(
            finite_hessian_gi, analytic_hessian_g[unsigned(i)]));
    }
}

// TEST_CASE("Displacment optimization test",
// "[!shouldfail][opt][displacements]")
// {
//     Eigen::MatrixX2d V;
//     V.resize(4, 2);
//     Eigen::MatrixX2d U;
//     U.resize(4, 2);
//     Eigen::MatrixX2i E;
//     E.resize(2, 2);
//     double volume_epsilon = 1.0;
//     DetectionMethod ccd_detection_method = DetectionMethod::BRUTE_FORCE;
//     Eigen::MatrixX2d Uopt;
//     Uopt.resize(4, 2);
//     Eigen::MatrixX2d actual_Uopt;
//     actual_Uopt.resize(4, 2);
//     OptimizationMethod opt_method;
//     unsigned max_iter = 1000;
//
//     SECTION("Horizontal edge falling")
//     {
//         // TODO: Fix this test
//         SECTION("MMA Optimizer") { opt_method = OptimizationMethod::MMA; }
//         SECTION("SLSQP Optimizer") { opt_method = OptimizationMethod::SLSQP;
//         } SECTION("IP Optimizer") { opt_method = OptimizationMethod::IP; } V
//         << -1, 1, 1, 1, -2, 0, 2, 0; U << 0, -2, 0, -2, 0, 0, 0, 0; E << 0,
//         1, 2, 3; displacements_optimization(V, U, E, volume_epsilon,
//             ccd_detection_method, opt_method, max_iter, Uopt);
//         actual_Uopt << 0, -1, 0, -1, 0, 0, 0, 0;
//         CHECK((Uopt - actual_Uopt).squaredNorm() ==
//         Approx(0.0).margin(1e-8));
//     }
//     // SECTION("Vertical edge falling")
//     // {
//     //     V << 0, 2, 0, 1, -1, 0, 1, 0;
//     //     U << 0, -3, 0, -3, 0, 0, 0, 0;
//     //     E << 0, 1, 2, 3;
//     //     displacements_nlopt_step(
//     //         V, U, E, volume_epsilon, ccd_detection_method, Uopt);
//     //     V(2, 0) *= 2;
//     //     V(3, 0) *= 2;
//     //     displacements_nlopt_step(
//     //         V, U, E, volume_epsilon, ccd_detection_method, actual_Uopt);
//     //     CHECK((Uopt - actual_Uopt).squaredNorm() ==
//     //     Approx(0.0).margin(1e-8));
//     // }
// }
