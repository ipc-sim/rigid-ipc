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
using namespace ipc;
using namespace ipc::rigid;

RigidBody
simple(Eigen::MatrixXd& vertices, Eigen::MatrixXi& edges, Pose<double> velocity)
{
    static int id = 0;
    int ndof = Pose<double>::dim_to_ndof(vertices.cols());
    return RigidBody(
        vertices, edges, Pose<double>::Zero(vertices.cols()), velocity,
        /*force=*/Pose<double>::Zero(vertices.cols()), /*density=*/1.0,
        /*is_dof_fixed=*/VectorXb::Zero(ndof), /*oriented=*/false,
        /*group=*/id++);
}

TEST_CASE("2D Rigid Body Transform", "[RB][RB-transform]")
{
    // Test vertices positions for given rb position
    Eigen::MatrixXd vertices_t0(4, 2);
    Eigen::MatrixXi edges(4, 2);
    Pose<double> velocity = Pose<double>::Zero(vertices_t0.cols());
    Pose<double> rb_step = Pose<double>::Zero(vertices_t0.cols());

    Eigen::MatrixXd vertices_step(4, 2), expected(4, 2);

    vertices_t0 << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_step.position << 0.5, 0.5;
        rb_step.rotation << 0.0;
        vertices_step = rb_step.position.transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_step.position << 0.0, 0.0;
        rb_step.rotation << 0.5 * igl::PI;
        vertices_step << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
    }

    SECTION("Translation and Rotation Case")
    {
        rb_step.position << 0.5, 0.5;
        rb_step.rotation << 0.5 * igl::PI;
        vertices_step << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        vertices_step += rb_step.position.transpose().replicate(4, 1);
    }
    expected = vertices_t0 + vertices_step;

    auto rb = simple(vertices_t0, edges, velocity);
    Pose<double> gamma_t1(
        rb.pose.position + rb_step.position,
        rb.pose.rotation + rb_step.rotation);
    Eigen::MatrixXd actual = rb.world_vertices<double>(gamma_t1);
    CHECK((expected - actual).squaredNorm() < 1E-6);
}

// TODO: Add 3D RB test
