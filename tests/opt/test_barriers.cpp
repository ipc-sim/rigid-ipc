#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <finitediff.hpp>
#include <barrier/barrier.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Test barriers and their derivatives", "[opt][barrier]")
{
    double s = 1e-4;

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

    Eigen::VectorXd x = Eigen::VectorXd::Random(1);
    // dis = std::uniform_real_distribution<double>(1e-8, s - 0.1);

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        x, [&phi](const Eigen::VectorXd& x) { return phi(x[0]); }, fgrad);
    Eigen::VectorXd grad = x.unaryExpr(phi_gradient);
    CHECK(fd::compare_gradient(fgrad, grad));

    fd::finite_gradient(x,
        [&phi_gradient](
            const Eigen::VectorXd& x) { return phi_gradient(x[0]); },
        fgrad);
    Eigen::VectorXd hessian = x.unaryExpr(phi_hessian);
    CHECK(fd::compare_gradient(fgrad, hessian));
}
