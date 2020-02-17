#include <catch2/catch.hpp>

#include <ccd/interval_root_finder.hpp>
#include <logger.hpp>

TEST_CASE("Simple interval arithmetic", "[ccd][interval]")
{
    ccd::Interval i(0, 1), j(4, 5);
    CHECK(i.lower() <= i.upper());
    CHECK(j.lower() <= j.upper());

    ccd::Interval result;

    result = i + j;
    CHECK(result.lower() <= result.upper());

    result = i - j;
    CHECK(result.lower() <= result.upper());

    result = i * j;
    CHECK(result.lower() <= result.upper());

    result = 1.0 / j;
    CHECK(result.lower() <= result.upper());

    result = i / j;
    CHECK(result.lower() <= result.upper());

    result = i + 10.0;
    CHECK(result.lower() == i.lower() + 10);
    CHECK(result.upper() == i.upper() + 10);
}

TEST_CASE("Cosine interval arithmetic", "[ccd][interval]")
{
    ccd::Interval r;

    double shift;
    SECTION("No shift") { shift = 0; }
    SECTION("2π shift") { shift = 2 * M_PI; }
    SECTION("-2π shift") { shift = -2 * M_PI; }
    SECTION("100π shift") { shift = 100 * M_PI; }
    SECTION("-100π shift") { shift = -100 * M_PI; }

    CAPTURE(shift);

    r = cos(ccd::Interval(-1, 7) + shift);
    CHECK(r.lower() == -1.0);
    CHECK(r.upper() == 1.0);

    r = cos(ccd::Interval(2, 4) + shift);
    CHECK(r.upper() < 0);

    r = cos(ccd::Interval(0, 1) + shift);
    CHECK(r.lower() > 0);

    r = cos(ccd::Interval(1, 2) + shift);
    CHECK(r.lower() < 0);
    CHECK(r.upper() > 0);
}

TEST_CASE("Sine interval arithmetic", "[ccd][interval]")
{
    ccd::Interval r;

    double shift = M_PI_2;
    SECTION("No shift") { shift += 0; }
    SECTION("2π shift") { shift += 2 * M_PI; }
    SECTION("-2π shift") { shift += -2 * M_PI; }
    SECTION("100π shift") { shift += 100 * M_PI; }
    SECTION("-100π shift") { shift += -100 * M_PI; }

    CAPTURE(shift);

    r = sin(ccd::Interval(-1, 7) + shift);
    CHECK(r.lower() == -1.0);
    CHECK(r.upper() == 1.0);

    r = sin(ccd::Interval(2, 4) + shift);
    CAPTURE(
        (ccd::Interval(2, 4) + shift).lower(),
        (ccd::Interval(2, 4) + shift).upper());
    CHECK(r.upper() < 0);

    r = sin(ccd::Interval(0, 1) + shift);
    CHECK(r.lower() > 0);

    r = sin(ccd::Interval(1, 2) + shift);
    CHECK(r.lower() < 0);
    CHECK(r.upper() > 0);
}

TEST_CASE("Root of simple function", "[ccd][interval]")
{
    double yshift = GENERATE(-1.1, -1.0 - 1e-11, -1.0, -0.5, 0.0, 0.5, 1.0);
    double a = 4, b = -4, c = 1 + yshift;
    auto f = [&](const ccd::Interval& x) { return a * x * x + b * x + c; };
    ccd::Interval sol;
    bool found_root =
        ccd::interval_root_finder(f, ccd::Interval(0, 1), sol, 1e-12);
    CHECK(found_root == (yshift >= -1.0 && yshift <= 0.0));
    if (found_root) {
        double actual_sol = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
        CHECK(median(sol) == Approx(actual_sol).margin(1e-12));
    }
}

TEST_CASE("Root of trig function", "[ccd][interval]")
{
    double yshift = GENERATE(0);
    auto f = [&](const ccd::Interval& x) { return cos(2 * M_PI * x) + yshift; };
    ccd::Interval sol;
    bool found_root =
        ccd::interval_root_finder(f, ccd::Interval(0, 1), sol, 1e-12);
    CHECK(found_root);
    if (found_root) {
        double actual_sol = 0.25;
        CHECK(median(sol) == Approx(actual_sol).margin(1e-12));
    }
}

TEST_CASE("Root of simple function with constraints", "[ccd][interval]")
{
    double yshift = GENERATE(-1.1, -1.0 - 1e-11, -1.0, -0.5, 0.0, 0.5, 1.0);
    double a = 4, b = -4, c = 1 + yshift;
    auto f = [&](const ccd::Interval& x) { return a * x * x + b * x + c; };
    auto constraint_predicate = [](const ccd::Interval& x) {
        return x.lower() > 0.5;
    };
    ccd::Interval sol;
    bool found_root = ccd::interval_root_finder(
        f, constraint_predicate, ccd::Interval(0, 1), sol, 1e-12);
    CHECK(found_root == (yshift >= -1.0 && yshift <= 0.0));
    if (found_root) {
        double actual_sol = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
        CHECK(median(sol) == Approx(actual_sol).margin(1e-12));
    }
}
