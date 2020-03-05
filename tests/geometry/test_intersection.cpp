#include <catch2/catch.hpp>

#include <ccd/interval.hpp>
#include <geometry/intersection.hpp>
#include <geometry/projection.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ccd::geometry;

TEST_CASE("Point-segment interval intersection", "[intersection]")
{
    typedef Eigen::VectorX3<ccd::Interval> VectorX3I;

    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorX3I p = VectorX3I::Zero(dim);
    p.y() = ccd::Interval(expected_distance);
    VectorX3I s0 = VectorX3I::Zero(dim);
    s0.x() = ccd::Interval(-10);
    VectorX3I s1 = VectorX3I::Zero(dim);
    s1.x() = ccd::Interval(10);

    bool intersecting = is_point_along_segment(p, s0, s1);
    CHECK(intersecting == abs(expected_distance) <= 1e-12);
    if (intersecting) {
        ccd::Interval alpha = project_point_to_line(p, s0, (s1 - s0).eval());
        CHECK(boost::numeric::median(alpha) == Approx(0.5).margin(1e-12));
    }
}

TEST_CASE("Segment-segment intersection", "[intersection]")
{
    // TODO
}

TEST_CASE("Point-triangle intersection", "[intersection]")
{
    // TODO
}
