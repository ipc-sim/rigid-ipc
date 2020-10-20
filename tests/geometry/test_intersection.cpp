#include <catch2/catch.hpp>

#include <ipc/friction/closest_point.hpp>

#include <interval/interval.hpp>
#include <geometry/intersection.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ccd::geometry;

TEST_CASE("Point-edge interval intersection", "[intersection]")
{
    typedef Eigen::VectorX3<ccd::Interval> VectorX3I;

    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorX3I p = VectorX3I::Zero(dim);
    p.y() = ccd::Interval(expected_distance);
    VectorX3I e0 = VectorX3I::Zero(dim);
    e0.x() = ccd::Interval(-10);
    VectorX3I e1 = VectorX3I::Zero(dim);
    e1.x() = ccd::Interval(10);

    bool intersecting = is_point_along_edge(p, e0, e1);
    CAPTURE(dim, expected_distance);
    CHECK(intersecting == (abs(expected_distance) <= 1e-12));
    if (intersecting) {
        ccd::Interval alpha = ipc::point_edge_closest_point(p, e0, e1);
        CHECK(boost::numeric::median(alpha) == Approx(0.5).margin(1e-12));
    }
}

TEST_CASE("Edge-edge intersection", "[intersection]")
{
    // TODO
}

TEST_CASE("Point-triangle intersection", "[intersection]")
{
    // TODO
}
