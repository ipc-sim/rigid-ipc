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
    Eigen::MatrixXd& vertices, Eigen::MatrixXi& edges, Pose<double> velocity)
{
    int ndof = Pose<double>::dim_to_ndof(vertices.cols());
    return RigidBody::from_points(vertices, edges,
        /*pose=*/Pose<double>(vertices.cols()), velocity,
        /*density=*/1.0,
        /*is_dof_fixed=*/Eigen::VectorXb::Zero(ndof),
        /*oriented=*/false);
}

} // namespace test_utils

TEST_CASE("Rigid Body System Transform", "[RB][RB-System][RB-System-transform]")
{
    using namespace test_utils;
    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(4, 2);
    Pose<double> velocity(vertices.cols());
    Pose<double> rb_displ_1(vertices.cols()), rb_displ_2(vertices.cols());

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_displ_1.position << 0.5, 0.5;
        rb_displ_2.position << 1.0, 1.0;

        // expected displacements
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2)
            = rb_displ_1.position.transpose().replicate(4, 1);
        expected.block(4, 0, 4, 2)
            = rb_displ_2.position.transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_displ_1.rotation << 0.5 * M_PI;
        rb_displ_2.rotation << M_PI;
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2) << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        expected.block(4, 0, 4, 2) << 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
            -1.0;
    }

    using namespace ccd::physics;

    std::vector<RigidBody> rbs;
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    RigidBodyAssembler assembler;
    assembler.init(rbs);

    std::vector<Pose<double>> poses(2);
    poses[0] = rb_displ_1 + rbs[0].pose;
    poses[1] = rb_displ_2 + rbs[1].pose;

    /// compute displacements between current and given positions
    /// TODO: update test to not need displacements
    Eigen::MatrixXd actual
        = assembler.world_vertices(poses) - assembler.world_vertices();
    CHECK((expected - actual).squaredNorm() < 1E-6);
}
