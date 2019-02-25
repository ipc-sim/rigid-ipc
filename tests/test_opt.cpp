#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Core>
#include <catch.hpp>
// #include <nlopt.hpp>

#include <autodiff/finitediff.hpp>
#include <ccd/optimization.hpp>

using namespace ccd::opt;

// TEST_CASE("Simple tests of NLOPT", "[nlopt]")
// {
//     auto f = [](Eigen::VectorXd x) { return x.squaredNorm(); };
//     auto df = [](Eigen::VectorXd x) { return 2 * x; };
//     nlopt::opt opt(nlopt::algorithm::NLOPT_LD_SLSQP, 10);
// }

TEST_CASE("Test ϕ and its derivatives", "[opt]")
{
    std::function<double(double, double)> phi;
    std::function<double(double, double)> phi_gradient;
    std::function<double(double, double)> phi_hessian;
    SECTION("Test phi_spline")
    {
        phi = phi_spline;
        phi_gradient = phi_spline_gradient;
        phi_hessian = phi_spline_hessian;
    }
    SECTION("Test phi_log")
    {
        phi = phi_log;
        phi_gradient = phi_log_gradient;
        phi_hessian = phi_log_hessian;
    }
    SECTION("Test phi_hookean")
    {
        phi = phi_hookean;
        phi_gradient = phi_hookean_gradient;
        phi_hessian = phi_hookean_hessian;
    }

    // Will be used to obtain a seed for the random number engine
    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    std::uniform_real_distribution<double> dis(1e-8, 10);
    double s = dis(gen);

    Eigen::VectorXd x(1);
    dis = std::uniform_real_distribution<double>(1e-8, s - 0.1);
    x[0] = dis(gen);

    Eigen::VectorXd fgrad(1);
    ccd::finite_gradient(
        x, [&phi, &s](const Eigen::VectorXd& x) { return phi(x[0], s); },
        fgrad);
    auto grad = phi_gradient(x[0], s);
    CHECK(fgrad[0] == Approx(grad));
    ccd::finite_gradient(
        x,
        [&phi_gradient, &s](
            const Eigen::VectorXd& x) { return phi_gradient(x[0], s); },
        fgrad, ccd::AccuracyOrder::EIGHTH);
    auto hessian = phi_hessian(x[0], s);
    CHECK(fgrad[0] == Approx(hessian));
}

TEST_CASE("Simple tests of Newton's Method", "[opt]")
{
    auto f = [](const Eigen::VectorXd& x) { return x.squaredNorm(); };
    auto gradient = [](const Eigen::VectorXd& x) { return 2 * x; };
    auto hessian = [](const Eigen::VectorXd& x) {
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
        x0, f, gradient, hessian, [](const Eigen::VectorXd&) { return true; });
    CHECK(min.squaredNorm() == Approx(0).margin(1e-10));
    CHECK(f(min) == Approx(0).margin(1e-10));
}

TEST_CASE("Simple tests of Newton's Method with inequlity constraints", "[opt]")
{
    auto f = [](const Eigen::VectorXd& x) { return x.squaredNorm(); };
    auto f_gradient = [](const Eigen::VectorXd& x) { return 2 * x; };
    auto f_hessian = [](const Eigen::VectorXd& x) {
        return 2 * Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    auto g = [](const Eigen::VectorXd& x) { return x[0] - 1; };
    auto g_gradient = [](Eigen::VectorXd x) {
        Eigen::VectorXd dg = Eigen::VectorXd::Zero(x.rows());
        dg[0] = 1;
        return dg;
    };
    auto g_hessian = [](Eigen::VectorXd x) {
        return Eigen::MatrixXd::Zero(x.rows(), x.rows());
    };

    double s = 1e-6;
    auto E = [&f, &g, &s](const Eigen::VectorXd& x) {
        return f(x) + phi_spline(g(x), s);
    };
    auto E_gradient = [&](const Eigen::VectorXd& x) {
        return f_gradient(x) + phi_spline_gradient(g(x), s) * g_gradient(x);
    };
    auto E_hessian = [&](const Eigen::VectorXd& x) {
        return f_hessian(x)
            + phi_spline_hessian(g(x), s) * g_gradient(x)
            * g_gradient(x).transpose()
            + phi_spline_gradient(g(x), s) * g_hessian(x);
    };

    Eigen::VectorXd x0(1);
    x0[0] = 5;
    Eigen::VectorXd min = newtons_method(x0, E, E_gradient, E_hessian,
        [&g](const Eigen::VectorXd& x) { return g(x) >= 0; });
    CHECK(min(0) == Approx(1.0));
}
