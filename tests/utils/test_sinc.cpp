#include <catch2/catch.hpp>

#include <iostream>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <autodiff/autodiff_types.hpp>
#include <utils/sinc.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("double sinc", "[sinc][double]")
{
    CHECK(sinc(0) == Approx(1));
    CHECK(sinc(1e-8) == Approx(1));
    CHECK(sinc(igl::PI / 2) == Approx(2 / igl::PI));
    CHECK(sinc(igl::PI) == Approx(0).margin(1e-16));
}

TEST_CASE("interval sinc", "[sinc][interval]")
{
    // All of these cases fall within the monotonic bound on sinc, so they are
    // a tight bound (ignoring rounding).
    Interval x, expected_y;

    SECTION("y=[1, 1]")
    {
        SECTION("[0, 0]") { x = Interval(0); };
        SECTION("[0, 1e-8]") { x = Interval(0, 1e-8); };
        SECTION("[-1e-8, 0]") { x = Interval(-1e-8, 0); };
        SECTION("[-1e-8, 1e-8]") { x = Interval(-1e-8, 1e-8); };
        expected_y = Interval(1);
    }
    SECTION("y=[0, 1]")
    {
        SECTION("[0, π]") { x = Interval(0, igl::PI); }
        SECTION("[-π, 0]") { x = Interval(-igl::PI, 0); }
        SECTION("[-π, π]") { x = Interval(-igl::PI, igl::PI); }
        expected_y = Interval(0, 1);
    }
    SECTION("y=[2/pi, 1]")
    {
        double PI_2 = igl::PI / 2;
        SECTION("[0, π/2]") { x = Interval(0, PI_2); }
        SECTION("[-π/2, 0]") { x = Interval(-PI_2, 0); }
        SECTION("[-π/2, π/2]") { x = Interval(-PI_2, PI_2); }
        expected_y = Interval(2 / igl::PI, 1);
    }

    CAPTURE(x.lower(), x.upper());
    Interval y = sinc(x);
    CHECK(y.lower() == Approx(expected_y.lower()).margin(1e-8));
    CHECK(y.upper() == Approx(expected_y.upper()).margin(1e-8));
}

TEST_CASE("interval sinc with looser bounds", "[sinc][interval]")
{
    Interval x, expected_y;

    SECTION("Non-monotonic")
    {
        double a = GENERATE(-5, -1, 0);
        x = Interval(a, 5);
        expected_y = Interval(-0.217233628211, 1);
    }
    SECTION("Monotonic outside bounds")
    {
        x = Interval(5, 7);
        expected_y = Interval(sinc(5), sinc(7));
    }
    SECTION("Monotonic far outside bounds")
    {
        x = Interval(21, 22);
        expected_y = Interval(sinc(22), sinc(21));
    }

    Interval y = sinc(x);
    CAPTURE(
        x.lower(), x.upper(), y.lower(), y.upper(), expected_y.lower(),
        expected_y.upper());
    CHECK(subset(expected_y, y));
    // Lossest bound
    CHECK(subset(y, Interval(-0.23, 1)));
    if (x.lower() > 0) {
        // Tighter bound if x.lower() > 1
        CHECK(subset(y, Interval(-1 / x.lower(), 1 / x.lower())));
    }
}

TEST_CASE("interval sinc_normx", "[sinc][interval]")
{
    VectorMax3<Interval> x = VectorMax3<Interval>::Zero(3);
    Interval expected_y;

    SECTION("Zero")
    {
        // x is already zero
        expected_y = Interval(1);
    }
    SECTION("Positive")
    {
        x(1) = Interval(igl::PI);
        expected_y = Interval(0);
    }
    SECTION("Negative")
    {
        x(1) = Interval(-igl::PI);
        expected_y = Interval(0);
    }
    SECTION("Mixed")
    {
        x(1) = Interval(-igl::PI, igl::PI);
        expected_y = Interval(0, 1);
    }

    Interval y = sinc_normx(x);
    CAPTURE(y.lower(), y.upper(), expected_y.lower(), expected_y.upper());
    CHECK(expected_y.lower() == Approx(y.lower()).margin(1e-8));
    CHECK(expected_y.upper() == Approx(y.upper()).margin(1e-8));
}

TEST_CASE("∇sinc(||x||)", "[sinc][vector][diff]")
{
    double sign = GENERATE(-1, 1);
    double val = GENERATE(0, 1e-8, igl::PI / 2, igl::PI, 5, 20, 100);
    int index = GENERATE(0, 1, 2);
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    x(index) = sign * val;

    Eigen::VectorXd fgrad(3);
    fd::finite_gradient(x, sinc_normx<double>, fgrad);

    Eigen::Vector3d grad = sinc_normx_grad(x);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("∇²sinc(||x||)", "[sinc][vector][diff]")
{
    double sign = GENERATE(-1, 1);
    double val = GENERATE(0, 1e-8, igl::PI / 2, igl::PI, 5, 20, 100);
    int index = GENERATE(0, 1, 2);
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    x(index) = sign * val;

    Eigen::MatrixXd fhess(3, 3);
    fd::finite_hessian(x, sinc_normx<double>, fhess);

    Eigen::Matrix3d hess = sinc_normx_hess(x);
    CHECK(fd::compare_hessian(hess, fhess));
}

template <typename T> VectorX<T> identity(const VectorX<T>& x) { return x; }
template <typename T> VectorX<T> square(const VectorX<T>& x)
{
    return x.array().square();
}
template <typename T> VectorX<T> const_sum(const VectorX<T>& x)
{
    return VectorX<T>::Constant(x.size(), x.sum());
}
template <typename T> VectorX<T> const_prod(const VectorX<T>& x)
{
    return VectorX<T>::Constant(x.size(), x.prod());
}

TEST_CASE("autodiff sinc(||x||)", "[sinc][vector][autodiff]")
{
    typedef AutodiffType<Eigen::Dynamic> Diff;
    int size = GENERATE(3, 6, 12);
    Diff::activate(size);

    double sign = GENERATE(-1, 1);
    double val = GENERATE(0, 1e-8, igl::PI / 2, igl::PI, 5, 20, 100);
    int index = GENERATE(0, 1, 2);
    Eigen::VectorXd x =
        GENERATE_COPY(Eigen::VectorXd::Zero(size), Eigen::VectorXd::Ones(size));
    x(size - 3 + index) = sign * val;

    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> g;
    std::function<Diff::D1VectorXd(const Diff::D1VectorXd&)> g_autodiff1;
    std::function<Diff::D2VectorXd(const Diff::D2VectorXd&)> g_autodiff2;
    SECTION("No chain rule")
    {
        g = identity<double>;
        g_autodiff1 = identity<Diff::DDouble1>;
        g_autodiff2 = identity<Diff::DDouble2>;
    }
    SECTION("x->x^2")
    {
        g = square<double>;
        g_autodiff1 = square<Diff::DDouble1>;
        g_autodiff2 = square<Diff::DDouble2>;
    }
    SECTION("x->const(sum(x))")
    {
        g = const_sum<double>;
        g_autodiff1 = const_sum<Diff::DDouble1>;
        g_autodiff2 = const_sum<Diff::DDouble2>;
    }
    SECTION("x->const(prod(x))")
    {
        g = const_prod<double>;
        g_autodiff1 = const_prod<Diff::DDouble1>;
        g_autodiff2 = const_prod<Diff::DDouble2>;
    }

    // Compute the expected values
    double expected_value = sinc_normx(VectorMax3d(g(x).tail<3>()));

    Eigen::VectorXd fgrad(size);
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& x) {
            return sinc_normx(VectorMax3d(g(x).tail<3>()));
        },
        fgrad);

    Eigen::MatrixXd fhess(size, size);
    fd::finite_hessian(
        x,
        [&](const Eigen::VectorXd& x) {
            return sinc_normx(VectorMax3d(g(x).tail<3>()));
        },
        fhess);

    // Check the DScalar1 version
    {
        Diff::D1VectorXd x_diff = Diff::d1vars(0, x);
        Diff::DDouble1 y = sinc_normx(
            VectorMax3<Diff::DDouble1>(g_autodiff1(x_diff).tail<3>()));

        double value = y.getValue();
        CHECK(value == Approx(expected_value));

        Eigen::VectorXd grad = y.getGradient();
        CHECK(fd::compare_gradient(grad, fgrad));
    }

    // Check the DScalar2 version
    {
        Diff::D2VectorXd x_diff = Diff::d2vars(0, x);
        Diff::DDouble2 y = sinc_normx(
            VectorMax3<Diff::DDouble2>(g_autodiff2(x_diff).tail<3>()));

        double value = y.getValue();
        CHECK(value == Approx(expected_value));

        Eigen::VectorXd grad = y.getGradient();
        CHECK(fd::compare_gradient(grad, fgrad));

        Eigen::MatrixXd hess = y.getHessian();
        CHECK(fd::compare_hessian(hess, fhess, 1e-2));
    }
}
