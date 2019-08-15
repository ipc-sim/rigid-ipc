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
    Eigen::MatrixXd actual = rb.world_vertices(gamma_t1);
    CHECK((expected - actual).squaredNorm() < 1E-6);
}

