#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <barrier/barrier.hpp>
#include <finitediff.hpp>

TEST_CASE("Test barriers and their derivatives", "[opt][barrier]")
{
    using namespace ccd::opt;
    double s = GENERATE(1e-4, 1, 10);

    BarrierType barrier_type =
        GENERATE(BarrierType::IPC, BarrierType::POLY_LOG, BarrierType::SPLINE);

    Eigen::VectorXd x = (Eigen::ArrayXd::Random(1) + 1) / 2 * s; // ∈ [0, s]

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& x) {
            return barrier(x[0], s, barrier_type);
        },
        fgrad);
    Eigen::VectorXd grad = x.unaryExpr([&](const double& xi) {
        return barrier_gradient(xi, s, barrier_type);
    });
    CAPTURE(barrier_type, s, x[0]);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Test spline barrier", "[opt][barrier]")
{
    using namespace ccd::opt;
    double s = GENERATE(1e-4, 1, 10);

    std::function<double(double)> phi;
    std::function<double(double)> phi_gradient;
    std::function<double(double)> phi_hessian;
    SECTION("Test spline_barrier")
    {
        phi = [&s](double x) { return spline_barrier(x, s); };
        phi_gradient = [&s](double x) { return spline_barrier_gradient(x, s); };
        phi_hessian = [&s](double x) { return spline_barrier_hessian(x, s); };
    }

    Eigen::VectorXd x(1);
    x << 0.5 * s; // ∈ [0, s]

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        x, [&phi](const Eigen::VectorXd& x) { return phi(x[0]); }, fgrad);
    Eigen::VectorXd grad = x.unaryExpr(phi_gradient);
    CAPTURE(s, x[0]);
    CHECK(fd::compare_gradient(fgrad, grad));

    fd::finite_gradient(
        x,
        [&phi_gradient](const Eigen::VectorXd& x) {
            return phi_gradient(x[0]);
        },
        fgrad);
    Eigen::VectorXd hessian = x.unaryExpr(phi_hessian);
    CHECK(fd::compare_gradient(fgrad, hessian));
}
