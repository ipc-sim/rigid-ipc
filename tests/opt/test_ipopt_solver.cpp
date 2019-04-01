#ifdef BUILD_WITH_IPOPT

#include <iomanip>
#include <iostream>

#include <catch.hpp>

#include <opt/ipopt_solver.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("IPOPT_hello_world", "[opt][Ipopt][Ipopt-Interface]")
{
    using namespace ccd::opt;

    // Our example problem has 4 variables (n), and 2 constraints (m).
    Eigen::VectorXd x0(4), x_lower(4), x_upper(4);
    Eigen::VectorXd g_lower(2), g_upper(2);

    x_lower.setConstant(1.0);
    x_upper.setConstant(5.0);

    int num_constraints = 2;
    g_lower[0] = 25;
    g_upper[0] = NO_UPPER_BOUND;
    g_lower[1] = g_upper[1] = 40.0;

    x0 << 1.0, 5.0, 5.0, 0.0;

    callback_f f = [](const Eigen::VectorXd& x) -> double {
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
    };

    callback_grad_f grad_f = [](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd grad_f(x.rows());
        grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        grad_f[1] = x[0] * x[3];
        grad_f[2] = x[0] * x[3] + 1;
        grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
        return grad_f;
    };

    callback_g g = [](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd g(2);
        g[0] = x[0] * x[1] * x[2] * x[3];
        g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        return g;
    };

    callback_jac_g jac_g = [](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd jac_g(2, 4);
        jac_g(0, 0) = x[1] * x[2] * x[3];
        jac_g(0, 1) = x[0] * x[2] * x[3];
        jac_g(0, 2) = x[0] * x[1] * x[3];
        jac_g(0, 3) = x[0] * x[1] * x[2];

        jac_g(1, 0) = 2 * x[0];
        jac_g(1, 1) = 2 * x[1];
        jac_g(1, 2) = 2 * x[2];
        jac_g(1, 3) = 2 * x[3];

        return jac_g;
    };

    Eigen::VectorXd expected(4);
    expected[0] = 1.0;
    expected[1] = 4.742999643576966129;
    expected[2] = 3.821149978948624959;
    expected[3] = 1.379408293215359738;

    OptimizationProblem problem(x0, f, grad_f, x_lower, x_upper,
        num_constraints, g, jac_g, g_lower, g_upper);
    SolverSettings settings;
    settings.verbosity = 0;

    auto result = ccd::opt::minimize_ipopt(problem, settings);
    CHECK((result.x - expected).norm() == Approx(0.0).margin(1e-8));
}

TEST_CASE("IPOPT_quadratic_linear_cnstr", "[opt][Ipopt][Ipopt-Interface]")
{
    using namespace ccd::opt;

    // Our example problem has 2 variables (n), and 0 constraints (m).
    // Example 16.4 from Nocedal, J, and S J Wright. 2006. Numerical
    // Optimization. Springer New York.
    int num_vars = 2, num_constraints = 3;

    Eigen::VectorXd x0(num_vars), x_lower(num_vars), x_upper(num_vars);
    Eigen::VectorXd g_lower(num_constraints), g_upper(num_constraints);

    x_lower.setConstant(0.0);
    x_upper.setConstant(NO_UPPER_BOUND);

    g_lower.setConstant(0.0);
    g_upper.setConstant(NO_UPPER_BOUND);

    x0 << 2.0, 0.0;

    callback_f f = [](const Eigen::VectorXd& x) -> double {
        // (x0 - 1)^2 + (x1 - 2.5)^2
        return (x[0] - 1) * (x[0] - 1) + (x[1] - 2.5) * (x[1] - 2.5);
    };

    callback_grad_f grad_f = [](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd grad_f(x.rows());
        grad_f[0] = 2 * (x[0] - 1);
        grad_f[1] = 2 * (x[1] - 2.5);
        return grad_f;
    };

    callback_g g
        = [num_constraints](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd g(num_constraints);
        g[0] = x[0] - 2 * x[1] + 2;
        g[1] = -x[0] - 2 * x[1] + 6;
        g[2] = -x[0] + 2 * x[1] + 2;
        return g;
    };

    callback_jac_g jac_g
        = [num_constraints](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd jac_g(num_constraints, x.rows());
        jac_g(0, 0) = 1.0;
        jac_g(0, 1) = -2.0;

        jac_g(1, 0) = -1.0;
        jac_g(1, 1) = -2.0;

        jac_g(2, 0) = -1.0;
        jac_g(2, 1) = 2.0;

        return jac_g;
    };

    Eigen::VectorXd expected(2);
    expected[0] = 1.4;
    expected[1] = 1.7;

    OptimizationProblem problem(x0, f, grad_f, x_lower, x_upper,
        num_constraints, g, jac_g, g_lower, g_upper);
    SolverSettings settings;
    settings.verbosity = 0;

    auto result = ccd::opt::minimize_ipopt(problem, settings);

    CHECK((result.x - expected).norm() < 1e-6);
}
TEST_CASE("IPOPT_quadratic_no_cnstr", "[opt][Ipopt][Ipopt-Interface]")
{
    using namespace ccd::opt;

    int num_vars = 2, num_constraints = 0;

    Eigen::VectorXd x0(num_vars), x_lower(num_vars), x_upper(num_vars);
    Eigen::VectorXd g_lower, g_upper;

    x_lower.setConstant(0.0);
    x_upper.setConstant(NO_UPPER_BOUND);

    x0 << 3.0, 2.0;

    callback_f f = [](const Eigen::VectorXd& x) -> double {
        // x1^2 + (x2 + 1)^2
        return (x[0] + 1) * (x[0] + 1) + (x[1] + 1) * (x[1] + 1);
    };

    callback_grad_f grad_f = [](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd grad_f(x.rows());
        grad_f[0] = 2 * (x[0] + 1);
        grad_f[1] = 2 * (x[1] + 1);
        return grad_f;
    };

    callback_g g
        = [num_constraints](const Eigen::VectorXd & /*x*/) -> Eigen::VectorXd {
        Eigen::VectorXd g(num_constraints);
        return g;
    };

    callback_jac_g jac_g
        = [num_constraints](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd jac_g(num_constraints, x.rows());
        return jac_g;
    };

    callback_intermediate callback
        = [](const Eigen::VectorXd& x, const double obj_value,
              const Eigen::VectorXd& dual, const double gamma, const int iteration) -> void {
        std::cout << "it=" << iteration << " x=" << x.transpose()
                  << " obj=" << obj_value << " y=" << dual.transpose()
                  << std::endl;
    };

    Eigen::VectorXd expected(2);
    expected[0] = 0.0;
    expected[1] = 0.0;

    OptimizationProblem problem(x0, f, grad_f, x_lower, x_upper,
        num_constraints, g, jac_g, g_lower, g_upper);
    SolverSettings settings;
    settings.verbosity = 0;
    settings.max_iter = 3000;

    auto result = ccd::opt::minimize_ipopt(problem, settings);
    CHECK((result.x - expected).norm() < 1e-6);
}
#endif
