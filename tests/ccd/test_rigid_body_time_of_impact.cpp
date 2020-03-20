#include <catch2/catch.hpp>

#include <ccd/rigid_body_time_of_impact.hpp>
#include <igl/PI.h>
#include <igl/edges.h>

#include <constants.hpp>

TEST_CASE("Rigid edge-vertex time of impact", "[ccd][rigid][rigid_toi]")
{
    using namespace ccd;
    using namespace ccd::physics;

    int dim = GENERATE(2);
    int ndof = Pose<double>::dim_to_ndof(dim);
    int pos_ndof = Pose<double>::dim_to_pos_ndof(dim);
    int rot_ndof = Pose<double>::dim_to_rot_ndof(dim);

    Eigen::MatrixXd bodyA_vertices(2, dim);
    bodyA_vertices.row(0) << -1, 0;
    bodyA_vertices.row(1) << 1, 0;
    Eigen::MatrixXd bodyB_vertices(2, dim);
    bodyB_vertices.row(0) << -2, 0;
    bodyB_vertices.row(1) << 2, 0;

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    int edge_id0 = GENERATE(0, 1); // Test both orders
    bodyB_edges.row(0) << edge_id0, edge_id0 == 0 ? 1 : 0;

    Pose<double> bodyA_pose(
        Eigen::Vector3d(0, 0.5, 0).head(pos_ndof),
        Eigen::Vector3d(0, 0, 0).head(rot_ndof));
    Pose<double> bodyB_pose = Pose<double>::Zero(dim);

    Pose<double> bodyA_velocity = Pose<double>::Zero(dim);
    Pose<double> bodyB_velocity = Pose<double>::Zero(dim);

    double expected_toi = -1;
    bool is_impact_expected = true;
    SECTION("Translation")
    {
        double posx = GENERATE(-2 - 1e-8, -2, -1, -0.5);
        double sign = GENERATE(-1, 1);
        bodyA_vertices.row(0).x() = sign * posx;
        bodyA_vertices.row(1).x() = -sign * posx;
        double vely = GENERATE(
            -1.5, -1.0 - 1e-12, -1.0, -1.0 + 1e-12, -0.5, -1e-12, 1e-12, 1.0);
        expected_toi = bodyA_pose.position.y() / -vely;
        bodyA_velocity.position = Eigen::Vector3d(0, vely, 0).head(pos_ndof);

        is_impact_expected =
            expected_toi >= 0 && expected_toi <= 1 && posx >= -2;
    }
    SECTION("Rotation")
    {
        double theta = igl::PI
            * GENERATE(-2, -7.0 / 6.0, -1, 0, 1.0 / 12.0, 1.0 / 6.0 - 1e-12,
                       1.0 / 6.0, 1.0 / 6.0 + 1e-12, 0.5, 1, 2, 100);
        expected_toi = (theta < 0 ? -7.0 : 1.0) / 6.0 * igl::PI / theta;
        is_impact_expected =
            theta >= igl::PI / 6.0 || theta <= -igl::PI * 7.0 / 6.0;
        bodyA_velocity.rotation = Eigen::Vector3d(0, 0, theta).tail(rot_ndof);
    }

    Eigen::VectorX6b bodyA_is_dof_fixed = Eigen::VectorX6b::Zero(ndof);
    Eigen::VectorX6b bodyB_is_dof_fixed = Eigen::VectorX6b::Ones(ndof);

    physics::RigidBody bodyA = physics::RigidBody::from_points(
        bodyA_vertices, bodyA_edges, bodyA_pose, bodyA_velocity, /*density=*/1,
        bodyA_is_dof_fixed, /*oriented=*/false);
    physics::RigidBody bodyB = physics::RigidBody::from_points(
        bodyB_vertices, bodyB_edges, bodyB_pose, bodyB_velocity, /*density=*/1,
        bodyB_is_dof_fixed, /*oriented=*/false);

    double toi;
    bool is_impacting = compute_edge_vertex_time_of_impact(
        bodyA, bodyA.pose, bodyA.velocity, /*vertex_id=*/0, bodyB, bodyB.pose,
        bodyB.velocity, /*edge_id=*/0, toi);
    CHECK(is_impacting == is_impact_expected);
    // if (is_impacting) {
    //     CHECK(
    //         toi
    //         == Approx(expected_toi)
    //                .margin(Constants::INTERVAL_ROOT_FINDER_TOL));
    // }
}

TEST_CASE("Rigid edge-edge time of impact", "[ccd][rigid][rigid_toi]") {}

TEST_CASE(
    "Rigid face-vertex time of impact", "[ccd][rigid][rigid_toi][face_vertex]")
{
    using namespace ccd;
    using namespace ccd::physics;
    int dim = 3;
    int ndof = Pose<double>::dim_to_ndof(dim);
    int pos_ndof = Pose<double>::dim_to_pos_ndof(dim);
    int rot_ndof = Pose<double>::dim_to_rot_ndof(dim);

    Eigen::MatrixXd bodyA_vertices(3, dim);
    bodyA_vertices.row(0) << -1, 0, 0;
    bodyA_vertices.row(1) << 1, 0, 0;
    bodyA_vertices.row(2) << 0, 1, 0;
    double y = GENERATE(1.0 + 1e-8, 1.0, 1.0 - 1e-8, 0.5, 1e-8, 0.0, 1e-8);
    Eigen::MatrixXd bodyB_vertices(3, dim);
    bodyB_vertices.row(0) << 0, y, 1;
    bodyB_vertices.row(1) << 1, y, 1.5;
    bodyB_vertices.row(2) << -1, y, 1.5;

    Eigen::MatrixXi bodyA_faces(1, 3);
    bodyA_faces.row(0) << 0, 1, 2;
    Eigen::MatrixXi bodyB_faces(1, 3);
    bodyB_faces.row(0) << 0, 1, 2;

    Eigen::MatrixXi bodyA_edges, bodyB_edges;
    igl::edges(bodyA_faces, bodyA_edges);
    igl::edges(bodyB_faces, bodyB_edges);

    Pose<double> bodyA_pose = Pose<double>::Zero(dim);
    Pose<double> bodyB_pose = Pose<double>::Zero(dim);

    Pose<double> bodyA_velocity = Pose<double>::Zero(dim);
    Pose<double> bodyB_velocity = Pose<double>::Zero(dim);

    double expected_toi = -1;
    bool is_impact_expected = false;
    double velz = GENERATE(-1.0 - 1e-8, -2.0, -10);
    SECTION("Translation")
    {
        expected_toi = -1 / velz;
        is_impact_expected = velz <= -1.0 && y >= 0 && y <= 1.0;
        bodyB_velocity.position.z() = velz;
    }
    // TODO: Rotation

    Eigen::VectorX6b bodyA_is_dof_fixed = Eigen::VectorX6b::Ones(ndof);
    Eigen::VectorX6b bodyB_is_dof_fixed = Eigen::VectorX6b::Zero(ndof);

    physics::RigidBody bodyA = physics::RigidBody::from_points(
        bodyA_vertices, bodyA_edges, bodyA_faces, bodyA_pose, bodyA_velocity,
        /*density=*/1, bodyA_is_dof_fixed, /*oriented=*/false);
    physics::RigidBody bodyB = physics::RigidBody::from_points(
        bodyB_vertices, bodyB_edges, bodyB_faces, bodyB_pose, bodyB_velocity,
        /*density=*/1, bodyB_is_dof_fixed, /*oriented=*/false);

    double toi;
    bool is_impacting = compute_face_vertex_time_of_impact(
        // Vertex body
        bodyB, bodyB.pose, bodyB.velocity, /*vertex_id=*/0,
        // Face body
        bodyA, bodyA.pose, bodyA.velocity, /*face_id=*/0,
        // Output time of impact
        toi);
    CAPTURE(y, velz);
    CHECK(is_impacting == is_impact_expected);
    if (is_impacting) {
        CHECK(
            toi
            == Approx(expected_toi)
                   .margin(Constants::INTERVAL_ROOT_FINDER_TOL));
    }
}
