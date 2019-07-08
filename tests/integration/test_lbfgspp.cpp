#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <solvers/lbfgspp/LBFGS.h>

#include <logger.hpp>

bool dummy_callback(const Eigen::VectorXd& x) { return true; }

TEST_CASE("Simple integration test of LBFGS++", "[opt][bfgs][lbfgspp][LBFGS++]")
{
    // f(x) = xᵀx / 2
    auto f = [](const Eigen::VectorXd& x, Eigen::VectorXd& grad) -> double {
        // ∇f(x) = x
        grad = x;
        return x.squaredNorm() / 2;
    };

    const int num_vars = 10;
    // Set up parameters
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-8;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSpp::LBFGSSolver<double> solver(param);

    // Initial guess
    Eigen::VectorXd x = Eigen::VectorXd::Random(num_vars);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(f, x, fx, dummy_callback);

    CHECK(x.squaredNorm() == Approx(0).margin(1e-8));
    CHECK(fx == Approx(0).margin(1e-8));
}

TEST_CASE(
    "Complex integration test of LBFGS++", "[opt][bfgs][lbfgspp][LBFGS++]")
{
    const int num_vars = 2;
    // f(x) = (a - x₀)² + b(x₁ - x₀²)²
    // https://en.wikipedia.org/wiki/Rosenbrock_function
    const double a = 1, b = 100;
    auto f = [a, b](const Eigen::VectorXd& x, Eigen::VectorXd& grad) -> double {
        assert(x.size() == num_vars);
        double tmp1 = (a - x(0));           // (a - x₀)
        double tmp2 = (x(1) - x(0) * x(0)); // (x₁ - x₀²)
        grad = Eigen::VectorXd::Zero(num_vars);
        // ∂f(x)/∂x₀ = -2(a - x₀) - 4b(x₁ - x₀²)x₀
        grad(0) = -2 * tmp1 - 4 * b * tmp2 * x(0);
        // ∂f(x)/∂x₁ = 2b(x₁ - x₀²)
        grad(1) = 2 * b * tmp2;
        return tmp1 * tmp1 + b * tmp2 * tmp2;
    };

    // Set up parameters
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-8;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSpp::LBFGSSolver<double> solver(param);

    // Initial guess
    Eigen::VectorXd x = Eigen::VectorXd::Random(num_vars);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(f, x, fx, dummy_callback);

    CHECK(x(0) == Approx(1).margin(1e-8));
    CHECK(x(1) == Approx(1).margin(1e-8));
}
