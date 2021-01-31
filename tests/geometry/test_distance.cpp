#include <catch2/catch.hpp>

#include <Eigen/Geometry>
#include <finitediff.hpp>
#include <igl/PI.h>

#include <ipc/distance/edge_edge.hpp>

#include <autodiff/autodiff_types.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ipc::rigid;

//-----------------------------------------------------------------------------
// Unsigned Distances
//-----------------------------------------------------------------------------

TEST_CASE("Edge-edge distance gradient", "[distance][gradient]")
{
    using namespace ipc::rigid;
    typedef AutodiffType<12> Diff;
    Diff::activate();

    // Generate a geometric space of
    double angle = 0;
    SECTION("Almost parallel")
    {
        double exponent = GENERATE(range(-6, 3));
        angle = pow(10, exponent) * igl::PI / 180.0;
    }
    // SECTION("Parallel") { angle = 0; }

    Diff::D2Vector3d ea0 = Diff::d2vars(0, Eigen::Vector3d(-1.0, 0, 0));
    Diff::D2Vector3d ea1 = Diff::d2vars(3, Eigen::Vector3d(+1.0, 0, 0));
    Diff::D2Vector3d eb0 =
        Diff::d2vars(6, Eigen::Vector3d(cos(angle), 1, sin(angle)));
    Diff::D2Vector3d eb1 = Diff::d2vars(
        9, Eigen::Vector3d(cos(angle + igl::PI), 1, sin(angle + igl::PI)));

    Diff::DDouble2 distance = ipc::edge_edge_distance(ea0, ea1, eb0, eb1);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(12);
    x.segment<3>(0) = Diff::get_value(ea0);
    x.segment<3>(3) = Diff::get_value(ea1);
    x.segment<3>(6) = Diff::get_value(eb0);
    x.segment<3>(9) = Diff::get_value(eb1);
    auto f = [](const Eigen::VectorXd& x) {
        return ipc::edge_edge_distance(
            Eigen::Vector3d(x.segment<3>(0)), Eigen::Vector3d(x.segment<3>(3)),
            Eigen::Vector3d(x.segment<3>(6)), Eigen::Vector3d(x.segment<3>(9)));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CAPTURE(angle, distance.getGradient().transpose(), fgrad.transpose());
    CHECK(distance.getValue() == Approx(1.0));
    CHECK(fd::compare_gradient(distance.getGradient(), fgrad));
    CHECK(distance.getHessian().squaredNorm() != 0.0);
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
    Eigen::Vector2d s0(-10, 0);
    Eigen::Vector2d s1(-9, 0);

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
