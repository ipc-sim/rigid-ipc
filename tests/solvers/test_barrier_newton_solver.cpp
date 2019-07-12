#include <catch2/catch.hpp>
#include <cmath>

#include <autodiff/finitediff.hpp>
#include <opt/ad_hoc_problem.hpp>
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

    double epsilon = GENERATE(1.0, 0.5, 1e-1, 5e-2);
    BarrierProblem barrier_problem(problem);

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
        Eigen::MatrixXd finite_hessian = barrier_problem.eval_hessian_f_approx(x);

        Eigen::MatrixXd analytic_hessian = barrier_problem.eval_hessian_f(x).toDense();
        CHECK(compare_jacobian(finite_hessian, analytic_hessian));

        CAPTURE(x, problem.g(x), epsilon,
            barrier_problem.eval_f(x), finite_grad, analytic_grad,
            finite_hessian, analytic_hessian);
    }
}


