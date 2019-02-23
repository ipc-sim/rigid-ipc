#include <Eigen/Core>
#include <catch.hpp>
// #include <nlopt.hpp>
#include <ccd/optimization.hpp>
#include <cmath>
#include <iostream>

using namespace ccd::opt;

// TEST_CASE("Simple tests of NLOPT", "[nlopt]")
// {
//     auto f = [](Eigen::VectorXd x) { return x.squaredNorm(); };
//     auto df = [](Eigen::VectorXd x) { return 2 * x; };
//     nlopt::opt opt(nlopt::algorithm::NLOPT_LD_SLSQP, 10);
// }

TEST_CASE("Simple tests of Newton's Method", "[opt]")
{
    auto f = [](Eigen::VectorXd x) { return x.squaredNorm(); };
    auto gradient = [](Eigen::VectorXd x) { return 2 * x; };
    auto hessian = [](Eigen::VectorXd x) {
        return 2 * Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };
    Eigen::VectorXd x0;
    SECTION("Dim = 1")
    {
        x0 = Eigen::VectorXd(1);
        x0 << 100;
    }
    SECTION("Dim = 3")
    {
        x0 = Eigen::VectorXd(3);
        x0 << 1, 2, 3;
    }
    SECTION("Dim = 10")
    {
        x0 = Eigen::VectorXd(10);
        x0 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
    }
    Eigen::VectorXd min = newtons_method(
        x0, f, gradient, hessian, [](Eigen::VectorXd) { return true; });
    CHECK(min.squaredNorm() == Approx(0).margin(1e-10));
    CHECK(f(min) == Approx(0).margin(1e-10));
}

TEST_CASE("Simple tests of Newton's Method with inequlity constraints", "[opt]")
{
    auto f = [](double x) { return x * x; };
    auto f_gradient = [](double x) { return 2 * x; };
    auto f_hessian = [](double x) { return 2; };

    auto g = [](double x) { return x - 1; };
    auto g_gradient = [](double x) { return 1; };
    auto g_hessian = [](double x) { return 0; };

    double si = 10;
    auto phi_hookean = [&g, &si](double x) { return pow(log(si * g(x)), 2); };
    auto phi_hookean_gradient = [&g, &g_gradient, &si](double x) {
        double gx = g(x);
        return 2 * log(si * gx) / gx * g_gradient(x);
    };
    auto phi_hookean_hessian = [&g, &g_gradient, &g_hessian, &si](double x) {
        double gx = g(x);
        double dgx = g_gradient(x);
        double f0 = 2 * log(si * gx);
        double df0 = 2 / gx * dgx;
        double f1 = gx;
        double df1 = dgx;
        double f2 = dgx;
        double df2 = g_hessian(x);
        return (f1 * (f2 * df0 + f0 * df2) - f0 * f2 * df1) / (f1 * f1);
    };

    auto E = [&f, &phi_hookean](
                 Eigen::VectorXd x) { return f(x(0)) + phi_hookean(x(0)); };
    auto E_gradient = [&f_gradient, &phi_hookean_gradient](Eigen::VectorXd x) {
        Eigen::VectorXd result(1);
        result << f_gradient(x(0)) + phi_hookean_gradient(x(0));
        return result;
    };
    auto E_hessian = [&](Eigen::VectorXd x) {
        Eigen::VectorXd result(1);
        result << f_hessian(x(0)) + phi_hookean_hessian(x(0));
        return result;
    };

    Eigen::VectorXd x0(1);
    // auto val = GENERATE(-10, -4, -1, 0, 1, 4, 10);
    // x0 << val;
    //
    // REQUIRE(f(x0) == Approx(val * val));
    // REQUIRE(f_gradient(x0)(0) == Approx(2 * val));
    // REQUIRE(f_hessian(x0).row(0)(0) == Approx(2));

    // REQUIRE(g(x0) == Approx(val));
    // REQUIRE(f_gradient(x0)(0) == Approx(2 * val));
    // REQUIRE(f_hessian(x0).row(0)(0) == Approx(2));

    x0(0) = 1.25;
    Eigen::VectorXd min = newtons_method(x0, E, E_gradient, E_hessian,
        [&g](Eigen::VectorXd x) { return g(x(0)) >= 0; });
    CHECK(min(0) == Approx(1.0906).epsilon(1e-4));
    // CHECK(f(min(0)) == Approx(0).margin(1e-6));
}
