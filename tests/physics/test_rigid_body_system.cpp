#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <physics/rigid_body_assembler.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

namespace test_utils {

using namespace ccd::physics;
RigidBody simple_rigid_body(
    Eigen::MatrixXd& vertices, Eigen::MatrixXi& edges, Eigen::MatrixXd velocity)
{
    return RigidBody::from_points(vertices, edges,
        /*mass=*/Eigen::VectorXd(),
        /*dof=*/Eigen::Vector3b::Zero(),
        /*oriented=*/false,
        /*position=*/Eigen::Vector3d::Zero(), velocity);
}

}

TEST_CASE("Rigid Body System Transform", "[RB][RB-System][RB-System-transform]")
{
    using namespace  test_utils;
    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(4, 2);
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d rb_displ_1, rb_displ_2;

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_displ_1 << 0.5, 0.5, 0.0;
        rb_displ_2 << 1.0, 1.0, 0.0;

        // expected displacements
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2)
            = rb_displ_1.segment(0, 2).transpose().replicate(4, 1);
        expected.block(4, 0, 4, 2)
            = rb_displ_2.segment(0, 2).transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_displ_1 << 0.0, 0.0, 0.5 * M_PI;
        rb_displ_2 << 0.0, 0.0, M_PI;
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2) << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        expected.block(4, 0, 4, 2) << 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
            -1.0;
    }

    using namespace ccd::physics;

    std::vector<RigidBody> rbs;
    RigidBodyAssembler assembler;
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    assembler.init(rbs);

    Eigen::VectorXd pos(6);
    pos.segment(0, 3) = rb_displ_1 + rbs[0].position;
    pos.segment(3, 3) = rb_displ_2 + rbs[1].position;

    /// compute displacements between current and given positions
    /// TODO: update test to not need displacements
    Eigen::MatrixXd actual = assembler.world_vertices(pos) - assembler.world_vertices();
    CHECK((expected - actual).squaredNorm() < 1E-6);

    //    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
    //    Eigen::DontAlignCols,
    //        ", ", ", ", "", "", " << ", ";");
    //    std::cout << expected.format(CommaInitFmt) << std::endl;
    //    std::cout << actual.format(CommaInitFmt) << std::endl;
}
