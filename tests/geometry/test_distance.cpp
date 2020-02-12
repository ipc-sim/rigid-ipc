#include <catch2/catch.hpp>

#include <geometry/distance.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ccd::geometry;

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

TEST_CASE("Point-point distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    Eigen::VectorX3d p0 = Eigen::VectorX3d::Zero(dim);
    Eigen::VectorX3d p1 = Eigen::VectorX3d::Zero(dim);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    SECTION("Aligned with X-axis") { p1(0) = expected_distance; }
    SECTION("Diagonal vector")
    {
        p1.setOnes();
        p1.normalize();
        p1 *= expected_distance;
    }
    double distance = point_point_distance<double>(p0, p1);
    CHECK(distance == Approx(abs(expected_distance)));
}

TEST_CASE("Point-segment distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = expected_distance;
    Eigen::VectorX3d s0 = Eigen::VectorX3d::Zero(dim);
    s0.x() = -10;
    Eigen::VectorX3d s1 = Eigen::VectorX3d::Zero(dim);
    s1.x() = 10;

    try {
        double distance = point_segment_distance<double>(p, s0, s1);
        CHECK(distance == Approx(abs(expected_distance)));
    } catch (ccd::NotImplementedError err) {
        // spdlog::error(err.what());
    }
}

TEST_CASE("Point-line distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = expected_distance;
    Eigen::VectorX3d s0 = Eigen::VectorX3d::Zero(dim);
    s0.x() = -9.9;
    Eigen::VectorX3d s1 = Eigen::VectorX3d::Zero(dim);
    s1.x() = -10;
    try {
        double distance = point_line_distance<double>(p, s0, s1);
        CHECK(distance == Approx(abs(expected_distance)));
    } catch (ccd::NotImplementedError err) {
        // spdlog::error(err.what());
    }
}

TEST_CASE("Segment-segment distance", "[distance]")
{
    // TODO
}

TEST_CASE("Point-triangle distance", "[distance]")
{
    // TODO
}

//-----------------------------------------------------------------------------
// Signed Distances
//-----------------------------------------------------------------------------

TEST_CASE("Point-line signed distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = expected_distance;
    Eigen::VectorX3d s0 = Eigen::VectorX3d::Zero(dim);
    s0.x() = -9.9;
    Eigen::VectorX3d s1 = Eigen::VectorX3d::Zero(dim);
    s1.x() = -10;
    try {
        double distance = point_line_signed_distance<double>(p, s0, s1);
        CHECK(distance == Approx(expected_distance));
    } catch (ccd::NotImplementedError err) {
        // spdlog::error(err.what());
    }
}

TEST_CASE("Line-line signed distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d line0_point0 = Eigen::VectorX3d::Zero(dim);
    line0_point0.x() = -9.9;
    line0_point0.y() = expected_distance;
    Eigen::VectorX3d line0_point1 = Eigen::VectorX3d::Zero(dim);
    line0_point0.x() = -10;
    line0_point1.y() = expected_distance;
    Eigen::VectorX3d line1_point0 = Eigen::VectorX3d::Zero(dim);
    line1_point0.x() = -9.9;
    Eigen::VectorX3d line1_point1 = Eigen::VectorX3d::Zero(dim);
    line1_point1.x() = -10;
    try {
        double distance = line_line_signed_distance<double>(
            line0_point0, line0_point1, line1_point0, line1_point1);
        CHECK(distance == Approx(abs(expected_distance)));
    } catch (ccd::NotImplementedError err) {
        // spdlog::error(err.what());
    }
}

TEST_CASE("Point-plane signed distance", "[distance]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = expected_distance;
    Eigen::VectorX3d p0 = Eigen::VectorX3d::Random(dim);
    p0.y() = 0;
    Eigen::VectorX3d n = Eigen::VectorX3d::Zero(dim);
    n.y() = 1;

    double distance = point_plane_signed_distance<double>(p, p0, n);
    CHECK(distance == Approx(expected_distance));
}
