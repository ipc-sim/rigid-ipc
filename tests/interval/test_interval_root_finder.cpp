#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <interval/interval_root_finder.hpp>
#include <constants.hpp>
#include <logger.hpp>

TEST_CASE("Root of simple function", "[ccd][interval]")
{
    double yshift = GENERATE(
        -1.1, -1.0 - 10 * ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL, -1.0,
        -0.5, -1e-4, 0.5, 1.0);
    double a = 4, b = -4, c = 1 + yshift;
    auto f = [&](const ipc::rigid::Interval& x) {
        return a * x * x + b * x + c;
    };
    ipc::rigid::Interval sol;
    bool found_root = ipc::rigid::interval_root_finder(
        f, ipc::rigid::Interval(0, 1),
        ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL, sol);
    CHECK(found_root == (yshift >= -1.0 && yshift <= 0.0));
    if (found_root) {
        double actual_sol = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
        CAPTURE(
            yshift, sol.lower(), sol.upper(), f(sol).lower(), f(sol).upper());
        CHECK(
            median(sol)
            == Approx(actual_sol)
                   .margin(ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL));
    }
}

TEST_CASE("Root of trig function", "[ccd][interval]")
{
    double yshift = GENERATE(0);
    auto f = [&](const ipc::rigid::Interval& x) {
        return cos(2 * igl::PI * x) + yshift;
    };
    ipc::rigid::Interval sol;
    bool found_root = ipc::rigid::interval_root_finder(
        f, ipc::rigid::Interval(0, 1),
        ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL, sol);
    CHECK(found_root);
    if (found_root) {
        double actual_sol = 0.25;
        CHECK(
            median(sol)
            == Approx(actual_sol)
                   .margin(ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL));
    }
}

TEST_CASE("Root of simple function with constraints", "[ccd][interval]")
{
    double yshift = GENERATE(-1.1, -1.0 - 1e-7, -1.0, -0.5, 0.0, 0.5, 1.0);
    double a = 4, b = -4, c = 1 + yshift;
    auto f = [&](const ipc::rigid::Interval& x) {
        return a * x * x + b * x + c;
    };
    auto constraint_predicate = [](const ipc::rigid::Interval& x) {
        return x.lower() > 0.5;
    };
    ipc::rigid::Interval sol;
    bool found_root = ipc::rigid::interval_root_finder(
        f, constraint_predicate, ipc::rigid::Interval(0, 1),
        ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL, sol,
        6 * ipc::rigid::Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);
    CHECK(found_root == (yshift >= -1.0 && yshift <= 0.0));
    if (found_root) {
        double actual_sol = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
        CHECK(
            median(sol)
            == Approx(actual_sol)
                   .margin(ipc::rigid::Constants::INTERVAL_ROOT_FINDER_TOL));
    }
}
