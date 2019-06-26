// Test the rigid body class.

#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <physics/rigid_body.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("Rigid Body Transform", "[RB][RB-transform]")
{
    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity;

    Eigen::MatrixX2d expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        velocity << 0.5, 0.5, 0.0;
        expected = velocity.segment(0, 2).transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        velocity << 0.0, 0.0, 0.5 * M_PI;
        expected << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
    }

    SECTION("Translation and Rotation Case")
    {
        velocity << 0.5, 0.5, 0.5 * M_PI;
        expected << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        expected += velocity.segment(0, 2).transpose().replicate(4, 1);
    }
    using namespace ccd::physics;

    auto rb = RigidBody::Centered(vertices, edges, velocity);
    Eigen::MatrixXd actual = rb.world_displacements();
    CHECK((expected - actual).squaredNorm() < 1E-6);
}

TEST_CASE("Rigid Body Gradient", "[RB][RB-gradient]")
{
    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity;

    Eigen::MatrixXd expected(8, 3);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case") { velocity << 0.5, 0.5, 0.0; }

    SECTION("90 Deg Rotation Case") { velocity << 0.0, 0.0, 0.5 * M_PI; }

    SECTION("Translation and Rotation Case")
    {
        velocity << 0.5, 0.5, 0.5 * M_PI;
    }

    for (int i = 0; i < 4; i++) {
        // x entry
        expected.row(i) << 1, 0,
            -vertices(i, 0) * sin(velocity.z())
            - vertices(i, 1) * cos(velocity.z());
        // y entry
        expected.row(4 + i) << 0, 1,
            vertices(i, 0) * cos(velocity.z())
            - vertices(i, 1) * sin(velocity.z());
    }

    using namespace ccd::physics;

    auto rb = RigidBody::Centered(vertices, edges, velocity);
    Eigen::MatrixXd actual = rb.world_displacements_gradient(rb.velocity);
    Eigen::MatrixXd finite
        = rb.world_displacements_gradient_finite(rb.velocity);
    Eigen::MatrixXd exact = rb.world_displacements_gradient_exact(rb.velocity);

    CHECK((expected - actual).squaredNorm() < 1E-6);
    CHECK((expected - finite).squaredNorm() < 1E-6);
    CHECK((expected - exact).squaredNorm() < 1E-6);
}

TEST_CASE("Rigid Body Hessian", "[RB][RB-hessian]")
{
    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity;

    std::array<Eigen::Matrix3d, 8> expected;

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case") { velocity << 0.5, 0.5, 0.0; }

    SECTION("90 Deg Rotation Case") { velocity << 0.0, 0.0, 0.5 * M_PI; }

    SECTION("Translation and Rotation Case")
    {
        velocity << 0.5, 0.5, 0.5 * M_PI;
    }

    for (int i = 0; i < 4; ++i) {
        Eigen::Matrix3d h = Eigen::Matrix3d::Zero(3, 3);
        h(2, 2) = -vertices(i, 0) * cos(velocity.z())
            + vertices(i, 1) * sin(velocity.z());
        expected[size_t(i)] = h;
    }
    for (int i = 0; i < 4; ++i) {
        Eigen::Matrix3d h = Eigen::Matrix3d::Zero(3, 3);
        h(2, 2) = -vertices(i, 0) * sin(velocity.z())
            - vertices(i, 1) * cos(velocity.z());
        expected[size_t(i) + 4] = h;
    }

    using namespace ccd::physics;

    auto rb = RigidBody::Centered(vertices, edges, velocity);
    std::vector<Eigen::Matrix3d> actual
        = rb.world_displacements_hessian(rb.velocity);
    std::vector<Eigen::Matrix3d> exact
        = rb.world_displacements_hessian_exact(rb.velocity);

    REQUIRE(expected.size() == actual.size());
    for (uint i = 0; i < 8; ++i) {
        CHECK((expected[i] - actual[i]).squaredNorm()
            == Approx(0.0).margin(1e-6));
        CHECK(
            (expected[i] - exact[i]).squaredNorm() == Approx(0.0).margin(1e-6));
    }
}
