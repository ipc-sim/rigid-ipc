#include <catch2/catch.hpp>
#include <cmath>

#include <autodiff/finitediff.hpp>
#include <solvers/barrier_solver.hpp>
#include <solvers/newton_solver.hpp>

#include <logger.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE(
    "Check barrier problem derivatives", "[opt][barrier][barrier_problem]")
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
    problem.hessian_g = [&problem](const Eigen::VectorXd& x)
        -> std::vector<Eigen::SparseMatrix<double>> {
        std::vector<Eigen::SparseMatrix<double>> hessian;
        for (long i = 0; i < problem.num_constraints; i++) {
            hessian.push_back(problem.hessian_f(x).sparseView());
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
    if (!std::isinf(barrier_problem.eval_f(x))) {
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
        finite_jacobian(x, barrier_problem.func_grad_f(), finite_hessian,
            AccuracyOrder::SECOND);
        Eigen::MatrixXd analytic_hessian = barrier_problem.eval_hessian_f(x);
        CHECK(compare_jacobian(finite_hessian, analytic_hessian));

        CAPTURE(x, problem.x_lower, problem.x_upper, problem.g(x),
            problem.g_lower, problem.g_upper, epsilon,
            barrier_problem.eval_f(x), finite_grad, analytic_grad,
            finite_hessian, analytic_hessian);
    }
}

TEST_CASE("Simple tests of Newton's Method with inequlity constraints",
    "[opt][barrier][barrier_solver]")
{
    spdlog::set_level(spdlog::level::off);
    // Setup solver
    BarrierSolver solver;
    solver.barrier_constraint = new BarrierConstraint();
    solver.barrier_constraint->barrier_epsilon = 1;
    solver.min_barrier_epsilon = 1e-12;
    solver.inner_solver.absolute_tolerance = 1e-12;
    // Setup problem
    int num_vars = 1, num_constraints = num_vars;
    AdHocProblem problem(num_vars, num_constraints);
    problem.is_dof_fixed = Eigen::MatrixXb::Zero(num_vars, 1);
    problem.x_lower = Eigen::VectorXd::Zero(num_vars);

    double x_shift = GENERATE(-1, 0, 1);

    problem.f = [x_shift](const Eigen::VectorXd& x) {
        return (x.array() + x_shift).matrix().squaredNorm() / 2;
    };
    problem.grad_f = [x_shift](const Eigen::VectorXd& x) {
        return (x.array() + x_shift).matrix();
    };
    problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    // TODO: Move the barrier function out of the constraint and back into the
    // barrier problem.
    problem.g = [](const Eigen::VectorXd& x) {
        // Eigen::VectorXd gx = Eigen::VectorXd::Zero(x.rows());
        // gx(0) = x(0);
        return Eigen::VectorXd::Zero(x.size());
    };
    problem.jac_g = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd dg = Eigen::MatrixXd::Zero(x.rows(), x.rows());
        // dg(0, 0) = 1;
        return dg;
    };
    problem.hessian_g = [](const Eigen::VectorXd& x) {
        return std::vector<Eigen::SparseMatrix<double>>(
            x.rows(), Eigen::SparseMatrix<double>(x.rows(), x.rows()));
    };

    problem.x0 = Eigen::VectorXd::Zero(num_vars);
    problem.x0(0) = GENERATE(0, 1, 5, 10);

    REQUIRE(problem.validate_problem());

    OptimizationResults results = solver.solve(problem);
    CAPTURE(x_shift, problem.x0(0), problem.f(problem.x0), results.success);
    if (problem.x0(0) + x_shift >= 0) {
        CHECK(results.x(0) == Approx(x_shift < 0 ? -x_shift : 0).margin(1e-6));
    }
    spdlog::set_level(spdlog::level::info);
}
