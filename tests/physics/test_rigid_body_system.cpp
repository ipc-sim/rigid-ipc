#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <physics/rigid_body_assembler.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------
using namespace ipc;
using namespace ipc::rigid;

RigidBody simple_rigid_body(
    Eigen::MatrixXd& vertices, Eigen::MatrixXi& edges, Pose<double> velocity)
{
    static int id = 0;
    int ndof = Pose<double>::dim_to_ndof(vertices.cols());
    return RigidBody(
        vertices, edges, /*pose=*/Pose<double>::Zero(vertices.cols()), velocity,
        /*force=*/Pose<double>::Zero(vertices.cols()), /*density=*/1.0,
        /*is_dof_fixed=*/VectorXb::Zero(ndof), /*oriented=*/false,
        /*group=*/id++);
}

TEST_CASE("Rigid Body System Transform", "[RB][RB-System][RB-System-transform]")
{
    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(4, 2);
    Pose<double> velocity = Pose<double>::Zero(vertices.cols());
    Pose<double> rb1_pose_t1 = Pose<double>::Zero(vertices.cols());
    Pose<double> rb2_pose_t1 = Pose<double>::Zero(vertices.cols());

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb1_pose_t1.position << 0.5, 0.5;
        rb2_pose_t1.position << 1.0, 1.0;

        // expected displacements
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2) =
            rb1_pose_t1.position.transpose().replicate(4, 1);
        expected.block(4, 0, 4, 2) =
            rb2_pose_t1.position.transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb1_pose_t1.rotation << 0.5 * igl::PI;
        rb2_pose_t1.rotation << igl::PI;
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2) << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        expected.block(4, 0, 4, 2) << 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
            -1.0;
    }

    std::vector<RigidBody> rbs;
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    rbs.push_back(simple_rigid_body(vertices, edges, velocity));
    RigidBodyAssembler assembler;
    assembler.init(rbs);

    Poses<double> poses = { { rb1_pose_t1, rb2_pose_t1 } };

    /// compute displacements between current and given positions
    /// TODO: update test to not need displacements
    Eigen::MatrixXd actual =
        assembler.world_vertices(poses) - assembler.world_vertices();
    CHECK((expected - actual).squaredNorm() < 1E-6);
}
