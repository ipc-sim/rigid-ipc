#include <catch2/catch.hpp>

#include <ccd/interval.hpp>
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
