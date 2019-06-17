#include <catch2/catch.hpp>

#include <solvers/newton_solver.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Simple tests of Newton's Method", "[opt][newtons method]")
{
    int num_vars = GENERATE(1, 10, 100), num_constraints = 0;
    // Setup solver
    NewtonSolver solver;
    solver.free_dof = Eigen::VectorXi::LinSpaced(num_vars, 0, num_vars - 1);
    // Setup problem
    // -----------------------------------------------------------------
    AdHocProblem problem(num_vars, num_constraints);
    problem.x0.setRandom();

    problem.f = [](const Eigen::VectorXd& x) { return x.squaredNorm() / 2.0; };
    problem.grad_f = [](const Eigen::VectorXd& x) { return x; };
    problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    REQUIRE(problem.validate_problem());

    OptimizationResults results = solver.solve(problem);
    REQUIRE(results.success);
    CHECK(results.x.squaredNorm() == Approx(0).margin(1e-6));
    CHECK(results.minf == Approx(0).margin(1e-6));
}
