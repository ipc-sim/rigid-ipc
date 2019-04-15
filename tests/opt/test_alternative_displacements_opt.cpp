#include <iostream>

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/alternative_displacement_opt.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Test the setup", "[opt][displacements]")
{
    Eigen::MatrixX2d V(4, 2);
    Eigen::MatrixX2d U(4, 2);
    Eigen::MatrixX2i E(2, 2);

    V << -1, 1, 1, 1, -2, 0, 2, 0;
    U << 0, -2, 0, -2, 0, 0, 0, 0;
    E << 0, 1, 2, 3;

    bool recompute_collision_set = GENERATE(false, true);
    DetectionMethod ccd_detection_method = DetectionMethod::BRUTE_FORCE;

    OptimizationProblem problem;
    alt::setup_displacement_optimization_problem(
        V, U, E, ccd_detection_method, recompute_collision_set, problem);

    CHECK(problem.num_vars == U.size());
    CHECK(problem.num_constraints == E.rows());

    CHECK(problem.x0.size() == U.size());
    CHECK(problem.x_lower.size() == U.size());
    CHECK(problem.x_upper.size() == U.size());
    CHECK(problem.g_lower.size() == E.rows());
    CHECK(problem.g_upper.size() == E.rows());

    Eigen::VectorXd x(problem.num_vars);
    x.setRandom();

    // Test ∇f
    Eigen::VectorXd finite_grad(problem.num_vars);
    finite_gradient(x, problem.f, finite_grad);
    Eigen::VectorXd analytic_grad = problem.grad_f(x);
    CHECK(compare_gradient(finite_grad, analytic_grad));

    // Test ∇²f
    Eigen::MatrixXd finite_hessian(problem.num_vars, problem.num_vars);
    finite_jacobian(x, problem.grad_f, finite_hessian);
    Eigen::MatrixXd analytic_hessian = problem.hessian_f(x);
    CHECK(compare_jacobian(finite_hessian, analytic_hessian));

    // Test ∇g
    Eigen::MatrixXd finite_jac(problem.num_constraints, problem.num_vars);
    finite_jacobian(x, problem.g, finite_jac);
    Eigen::MatrixXd analytic_jac = problem.jac_g(x);
    CHECK(compare_jacobian(finite_jac, analytic_jac));

    // Test ∇²g
    Eigen::MatrixXd finite_hessian_gi(problem.num_vars, problem.num_vars);
    finite_hessian_gi.setOnes();
    long i;
    auto diff_i = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        assert(i >= 0 && i < problem.num_constraints);
        auto jac = problem.jac_g(x);
        return jac.row(i);
    };
    std::vector<Eigen::MatrixXd> analytic_hessian_g = problem.hessian_g(x);
    for (i = 0; i < problem.num_constraints; i++) {
        finite_jacobian(x, diff_i, finite_hessian_gi);
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
