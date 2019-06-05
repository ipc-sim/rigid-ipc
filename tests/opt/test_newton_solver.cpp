#include <catch2/catch.hpp>

#include <opt/newtons_method.hpp>

using namespace ccd;
using namespace opt;

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

    Eigen::VectorXi free_dof(num_vars);
    for (int i = 0; i < free_dof.size(); i++) {
        free_dof(i) = i;
    }

    OptimizationResults results
        = newtons_method(problem, free_dof, 1e-6, 1e-6, 100);
    REQUIRE(results.success);
    CHECK(results.x.squaredNorm() == Approx(0).margin(1e-6));
    CHECK(results.minf == Approx(0).margin(1e-6));
}
