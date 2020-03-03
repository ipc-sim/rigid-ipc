#include <catch2/catch.hpp>

#include <geometry/distance.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ccd::geometry;

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

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

template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }

TEST_CASE("Point-line signed distance", "[distance]")
{
    double expected_distance = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector2d p = Eigen::Vector2d::Random();
    p.y() = expected_distance;
    Eigen::Vector2d s0(-9, 0);
    Eigen::Vector2d s1(-10, 0);

    double distance = point_line_signed_distance(p, s0, s1);
    CAPTURE(distance, expected_distance);
    CHECK(sign(distance) == sign(expected_distance));
}

TEST_CASE("Line-line signed distance", "[distance]")
{
    double expected_distance = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d line0_point0(-9.9, expected_distance, 0);
    Eigen::Vector3d line0_point1(-10, expected_distance, 0);
    Eigen::Vector3d line1_point0(0, 0, -10);
    Eigen::Vector3d line1_point1(0, 0, -9.9);

    double distance = line_line_signed_distance(
        line0_point0, line0_point1, line1_point0, line1_point1);
    CAPTURE(distance, expected_distance);
    CHECK(sign(distance) == sign(expected_distance));
}

TEST_CASE("Point-plane signed distance", "[distance]")
{
    double expected_distance = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d p = Eigen::Vector3d::Random();
    p.y() = expected_distance;
    Eigen::Vector3d p0 = Eigen::Vector3d::Random();
    Eigen::Vector3d p1 = Eigen::Vector3d::Random();
    Eigen::Vector3d p2 = Eigen::Vector3d::Random();
    p0.y() = p1.y() = p2.y() = 0;
    if (Eigen::Vector3d::UnitY().dot((p1 - p0).cross(p2 - p0)) < 0) {
        std::swap(p1, p2);
    }

    double distance = point_plane_signed_distance(p, p0, p1, p2);
    CAPTURE(distance, expected_distance);
    CHECK(sign(distance) == sign(expected_distance));
}
