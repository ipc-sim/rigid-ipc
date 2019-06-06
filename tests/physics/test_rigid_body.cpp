#include <iomanip>
#include <iostream>

#include <catch.hpp>

#include <igl/PI.h>

#include <physics/rigid_body_system.hpp>

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


    //    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
    //    Eigen::DontAlignCols,
    //        ", ", ", ", "", "", " << ", ";");
    //    std::cout << expected.format(CommaInitFmt) << std::endl;
    //    std::cout << actual.format(CommaInitFmt) << std::endl;
}
