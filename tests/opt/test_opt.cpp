#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Core>
#include <catch.hpp>
#include <nlopt.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

using namespace ccd::opt;

double test_func(const std::vector<double>& x, std::vector<double>& grad,
    void* /* my_func_data */)
{
    Eigen::VectorXd X
        = Eigen::Map<const Eigen::VectorXd>(x.data(), long(x.size()));
    if (!grad.empty()) {
        Eigen::VectorXd gradX = 2 * X;
        assert(size_t(gradX.size()) == size_t(grad.size()));
        for (int i = 0; i < int(grad.size()); i++) {
            grad[size_t(i)] = gradX[i];
        }
    }
    return X.squaredNorm();
}

TEST_CASE("Simple tests of NLOPT", "[opt][nlopt]")
{
    const int DIM = 10;
    nlopt::opt opt(nlopt::LD_MMA, DIM);
    std::vector<double> lb(DIM, -HUGE_VAL); // lower bounds
    lb[0] = 1;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(test_func, nullptr);
    opt.set_xtol_rel(1e-12);
    std::vector<double> x(DIM, 100); // some initial guess
    double minf;                     // the minimum objective value, upon return
    nlopt::result result = opt.optimize(x, minf);
    CHECK(x[0] == Approx(1.0));
    CHECK(minf == Approx(1.0));
}

TEST_CASE("Test ϕ and its derivatives", "[opt][barrier]")
{
    std::function<double(double, double)> phi;
    std::function<double(double, double)> phi_gradient;
    std::function<double(double, double)> phi_hessian;
    SECTION("Test spline_barrier")
    {
        phi = spline_barrier;
        phi_gradient = spline_barrier_gradient;
        phi_hessian = spline_barrier_hessian;
    }
    SECTION("Test log_barrier")
    {
        phi = log_barrier;
        phi_gradient = log_barrier_gradient;
        phi_hessian = log_barrier_hessian;
    }
    SECTION("Test hookean_barrier")
    {
        phi = hookean_barrier;
        phi_gradient = hookean_barrier_gradient;
        phi_hessian = hookean_barrier_hessian;
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
    Eigen::VectorXd x;
    SECTION("Dim = 1")
    {
        x = Eigen::VectorXd(1);
        x << 100;
    }
    SECTION("Dim = 3")
    {
        x = Eigen::VectorXd(3);
        x << 1, 2, 3;
    }
    SECTION("Dim = 10")
    {
        x = Eigen::VectorXd(10);
        x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
    }
    auto no_constraint = [](const Eigen::VectorXd&) { return true; };
    newtons_method(x, f, gradient, hessian, no_constraint);
    CHECK(x.squaredNorm() == Approx(0).margin(1e-10));
    CHECK(f(x) == Approx(0).margin(1e-10));
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
        return f(x) + spline_barrier(g(x), s);
    };
    auto E_gradient = [&](const Eigen::VectorXd& x) {
        return f_gradient(x) + spline_barrier_gradient(g(x), s) * g_gradient(x);
    };
    auto E_hessian = [&](const Eigen::VectorXd& x) {
        return f_hessian(x)
            + spline_barrier_hessian(g(x), s) * g_gradient(x)
            * g_gradient(x).transpose()
            + spline_barrier_gradient(g(x), s) * g_hessian(x);
    };

    Eigen::VectorXd x(1);
    x[0] = 5;
    newtons_method(x, E, E_gradient, E_hessian,
        [&g](const Eigen::VectorXd& x) { return g(x) >= 0; });
    CHECK(x[0] == Approx(1.0));
}
