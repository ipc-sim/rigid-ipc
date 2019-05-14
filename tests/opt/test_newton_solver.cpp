#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Core>
#include <catch.hpp>
#include <nlopt.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/barrier.hpp>
#include <opt/barrier_newton_solver.hpp>
#include <opt/newtons_method.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Test ϕ and its derivatives", "[opt][barrier]")
{
    // Will be used to obtain a seed for the random number engine
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    std::uniform_real_distribution<double> dis(1e-8, 10);
    double s = dis(gen);

    std::function<double(double)> phi;
    std::function<double(double)> phi_gradient;
    std::function<double(double)> phi_hessian;
    SECTION("Test spline_barrier")
    {
        phi = [&s](double x) { return spline_barrier(x, s); };
        phi_gradient = [&s](double x) { return spline_barrier_gradient(x, s); };
        phi_hessian = [&s](double x) { return spline_barrier_hessian(x, s); };
    }
    SECTION("Test log_barrier")
    {
        phi = [&s](double x) { return log_barrier(x, s); };
        phi_gradient = [&s](double x) { return log_barrier_gradient(x, s); };
        phi_hessian = [&s](double x) { return log_barrier_hessian(x, s); };
    }
    SECTION("Test hookean_barrier")
    {
        phi = [&s](double x) { return hookean_barrier(x, s); };
        phi_gradient
            = [&s](double x) { return hookean_barrier_gradient(x, s); };
        phi_hessian = [&s](double x) { return hookean_barrier_hessian(x, s); };
    }

    Eigen::VectorXd x(1);
    dis = std::uniform_real_distribution<double>(1e-8, s - 0.1);
    x[0] = dis(gen);

    Eigen::VectorXd fgrad(1);
    ccd::finite_gradient(
        x, [&phi](const Eigen::VectorXd& x) { return phi(x[0]); }, fgrad);
    Eigen::VectorXd grad = x.unaryExpr(phi_gradient);
    CHECK(compare_gradient(fgrad, grad));

    ccd::finite_gradient(
        x,
        [&phi_gradient](
            const Eigen::VectorXd& x) { return phi_gradient(x[0]); },
        fgrad);
    Eigen::VectorXd hessian = x.unaryExpr(phi_hessian);
    CHECK(compare_gradient(fgrad, hessian));
}

TEST_CASE("Simple tests of Newton's Method", "[opt][newtons method]")
{
    // Setup problem
    // -----------------------------------------------------------------
    int num_vars = GENERATE(1, 10, 100), num_constraints = 0;
    AdHocProblem problem(num_vars, num_constraints);
    problem.x0.setRandom();

    problem.f = [](const Eigen::VectorXd& x) { return x.squaredNorm() / 2.0; };
    problem.grad_f = [](const Eigen::VectorXd& x) { return x; };
    problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    REQUIRE(problem.validate_problem());

    OptimizationResults results = newtons_method(problem, 1e-6, 1e-6, 100);
    REQUIRE(results.success);
    CHECK(results.x.squaredNorm()
        == Approx(0).margin(1e-6));
    CHECK(results.minf == Approx(0).margin(1e-6));
}

TEST_CASE("Check barrier problem derivatives")
{
    // TODO: Generate a random problem and test the derivatives
    int num_vars = GENERATE(1, 2, 5, 10), num_constraints = 1;
    AdHocProblem problem(num_vars, num_constraints);

    problem.f = [](const Eigen::VectorXd& x) { return x.squaredNorm() / 2; };
    problem.grad_f = [](const Eigen::VectorXd& x) { return x; };
    problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    problem.g = [&problem](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd gx(problem.num_constraints);
        gx.setConstant(problem.f(x));
        return gx;
    };
    problem.jac_g = [&problem](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        return problem.grad_f(x).transpose();
    };
    problem.hessian_g = [&problem](const Eigen::VectorXd& x) {
        std::vector<Eigen::MatrixXd> hessian;
        for (long i = 0; i < problem.num_constraints; i++) {
            hessian.push_back(problem.hessian_f(x));
        }
        return hessian;
    };

    problem.x_lower.setConstant(GENERATE(double(NO_LOWER_BOUND), -2.0));
    problem.x_upper.setConstant(GENERATE(double(NO_UPPER_BOUND), 2.0));
    problem.g_lower.setConstant(GENERATE(double(NO_LOWER_BOUND), 0.5));
    problem.g_upper.setConstant(GENERATE(double(NO_UPPER_BOUND), 4.0));

    REQUIRE(problem.validate_problem());


    double epsilon = GENERATE(1.0, 0.5, 1e-1, 5e-2);
    BarrierProblem barrier_problem(problem, epsilon);
    REQUIRE(barrier_problem.validate_problem());

    Eigen::VectorXd x(num_vars);
    x.setConstant(GENERATE(-1.0, 1.0) * GENERATE(1e-1, 1.0, 2.0 - 1e-3, 4.0));

    // If the function evaluates to infinity then the finite differences will
    // not work. I assume in the definition of the barrier gradient that
    // d/dx ∞ = 0.
    if (!isinf(barrier_problem.eval_f(x))) {
        // Use a higher order finite difference method because the function near
        // the boundary becomes very non-linear. This problem worsens as the ϵ
        // of the boundary gets smaller.

        // Test ∇f
        Eigen::VectorXd finite_grad(barrier_problem.num_vars);
        finite_gradient(
            x, barrier_problem.func_f(), finite_grad, AccuracyOrder::SECOND);
        Eigen::VectorXd analytic_grad = barrier_problem.eval_grad_f(x);
        CHECK(compare_gradient(finite_grad, analytic_grad));

        // Test ∇²f
        Eigen::MatrixXd finite_hessian(
            barrier_problem.num_vars, barrier_problem.num_vars);
        finite_jacobian(
            x, barrier_problem.func_grad_f(), finite_hessian, AccuracyOrder::SECOND);
        Eigen::MatrixXd analytic_hessian = barrier_problem.eval_hessian_f(x);
        CHECK(compare_jacobian(finite_hessian, analytic_hessian));

        CAPTURE(x, problem.x_lower, problem.x_upper, problem.g(x),
            problem.g_lower, problem.g_upper, epsilon, barrier_problem.eval_f(x),
            finite_grad, analytic_grad, finite_hessian, analytic_hessian);
    }
}

TEST_CASE("Simple tests of Newton's Method with inequlity constraints",
    "[opt][barrier]")
{
    // Setup problem
    // -----------------------------------------------------------------
    int num_vars = 1, num_constraints = num_vars;
    AdHocProblem constrained_problem(num_vars, num_constraints);


    SECTION("Constraint is in g(x)") { constrained_problem.g_lower(0) = 1; }
    SECTION("Constraint is in x_lower") { constrained_problem.x_lower(0) = 1; }

    constrained_problem.f
        = [](const Eigen::VectorXd& x) { return x.squaredNorm() / 2; };
    constrained_problem.grad_f = [](const Eigen::VectorXd& x) { return x; };
    constrained_problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    constrained_problem.g = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd gx = Eigen::VectorXd::Zero(x.rows());
        gx(0) = x(0);
        return gx;
    };
    constrained_problem.jac_g = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd dg = Eigen::MatrixXd::Zero(x.rows(), x.rows());
        dg(0, 0) = 1;
        return dg;
    };
    constrained_problem.hessian_g = [](const Eigen::VectorXd& x) {
        std::vector<Eigen::MatrixXd> hessian;
        for (long i = 0; i < x.rows(); i++) {
            hessian.push_back(Eigen::MatrixXd::Zero(x.rows(), x.rows()));
        }
        return hessian;
    };

    REQUIRE(constrained_problem.validate_problem());

    double s = 1e-6;
    BarrierProblem unconstrained_problem(constrained_problem, s);

    unconstrained_problem.x0(0) = 5;

    REQUIRE(unconstrained_problem.validate_problem());


    OptimizationResults results
        = newtons_method(unconstrained_problem, 1e-6, 1e-6, 100);
    // REQUIRE(results.success);
    CHECK(results.x(0) == Approx(1.0).margin(1e-6));
}
