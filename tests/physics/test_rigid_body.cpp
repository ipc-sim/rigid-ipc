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
    // Test vertices positions for given rb position

    Eigen::MatrixX2d vertices_t0(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d rb_step;

    Eigen::MatrixX2d vertices_step(4, 2), expected(4, 2);

    vertices_t0 << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_step << 0.5, 0.5, 0.0;
        vertices_step = rb_step.segment(0, 2).transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_step << 0.0, 0.0, 0.5 * M_PI;
        vertices_step << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
    }

    SECTION("Translation and Rotation Case")
    {
        rb_step << 0.5, 0.5, 0.5 * M_PI;
        vertices_step << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        vertices_step += rb_step.segment(0, 2).transpose().replicate(4, 1);
    }
    using namespace ccd::physics;
    expected = vertices_t0 + vertices_step;

    auto rb = RigidBody::from_velocity(vertices_t0, edges, velocity);
    Eigen::VectorXd gamma_t1 = rb.position + rb_step;
    Eigen::MatrixXd actual = rb.world_vertices<double>(gamma_t1);
    CHECK((expected - actual).squaredNorm() < 1E-6);
}

TEST_CASE("Rigid Body Gradient", "[RB][RB-gradient]")
{
    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d rb_velocity = Eigen::Vector3d::Zero(), rb_displacement;

    Eigen::MatrixXd expected(8, 3);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case") { rb_displacement << 0.5, 0.5, 0.0; }

    SECTION("90 Deg Rotation Case") { rb_displacement << 0.0, 0.0, 0.5 * M_PI; }

    SECTION("Translation and Rotation Case")
    {
        rb_displacement << 0.5, 0.5, 0.5 * M_PI;
    }

    for (int i = 0; i < 4; i++) {
        // x entry
        expected.row(i) << 1, 0,
            -vertices(i, 0) * sin(rb_displacement.z())
            - vertices(i, 1) * cos(rb_displacement.z());
        // y entry
        expected.row(4 + i) << 0, 1,
            vertices(i, 0) * cos(rb_displacement.z())
            - vertices(i, 1) * sin(rb_displacement.z());
    }

    using namespace ccd::physics;
    auto rb = RigidBody::from_velocity(vertices, edges, rb_velocity);

    Eigen::VectorXd gamma_t1 = rb.position + rb_displacement;
    Eigen::MatrixXd actual = rb.world_vertices_gradient(gamma_t1);
    Eigen::MatrixXd finite = rb.world_vertices_gradient_finite(gamma_t1);
    //    Eigen::MatrixXd exact = rb.world_vertices_gradient_exact(q1);

    CHECK((expected - actual).squaredNorm() < 1E-6);
    CHECK((expected - finite).squaredNorm() < 1E-6);
    //    CHECK((expected - exact).squaredNorm() < 1E-6);
}

TEST_CASE("Rigid Body Hessian", "[RB][RB-hessian]")
{
    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d rb_velocity = Eigen::Vector3d::Zero(), rb_displacement;

    std::array<Eigen::Matrix3d, 8> expected;

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case") { rb_displacement << 0.5, 0.5, 0.0; }

    SECTION("90 Deg Rotation Case") { rb_displacement << 0.0, 0.0, 0.5 * M_PI; }

    SECTION("Translation and Rotation Case")
    {
        rb_displacement << 0.5, 0.5, 0.5 * M_PI;
    }

    for (int i = 0; i < 4; ++i) {
        Eigen::Matrix3d h = Eigen::Matrix3d::Zero(3, 3);
        h(2, 2) = -vertices(i, 0) * cos(rb_displacement.z())
            + vertices(i, 1) * sin(rb_displacement.z());
        expected[size_t(i)] = h;
    }
    for (int i = 0; i < 4; ++i) {
        Eigen::Matrix3d h = Eigen::Matrix3d::Zero(3, 3);
        h(2, 2) = -vertices(i, 0) * sin(rb_displacement.z())
            - vertices(i, 1) * cos(rb_displacement.z());
        expected[size_t(i) + 4] = h;
    }

    using namespace ccd::physics;
    auto rb = RigidBody::from_velocity(vertices, edges, rb_velocity);

    Eigen::VectorXd gamma_t1 = rb.position + rb_displacement;
    std::vector<Eigen::Matrix3d> actual = rb.world_vertices_hessian(gamma_t1);
    std::vector<Eigen::Matrix3d> exact
        = rb.world_vertices_hessian_exact(gamma_t1);

    REQUIRE(expected.size() == actual.size());
    for (uint i = 0; i < 8; ++i) {
        CHECK((expected[i] - actual[i]).squaredNorm()
            == Approx(0.0).margin(1e-6));
        CHECK(
            (expected[i] - exact[i]).squaredNorm() == Approx(0.0).margin(1e-6));
    }
}
