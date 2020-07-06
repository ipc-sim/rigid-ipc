#include <catch2/catch.hpp>

#include <ccd/rigid_body_time_of_impact.hpp>
#include <igl/PI.h>
#include <igl/edges.h>

#include <constants.hpp>

using namespace ccd;
using namespace ccd::physics;

RigidBody create_body(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
{
    static int id = 0;
    int dim = vertices.cols();
    Pose<double> pose = Pose<double>::Zero(dim);
    RigidBody rb = RigidBody::from_points(
        vertices, edges, faces, pose,
        /*velocity=*/Pose<double>::Zero(pose.dim()),
        /*force=*/Pose<double>::Zero(pose.dim()),
        /*denisty=*/1.0,
        /*is_dof_fixed=*/Eigen::VectorX6b::Zero(pose.ndof()),
        /*oriented=*/false,
        /*group_id=*/id++);
    rb.vertices = vertices; // Cancel out the inertial rotation for testing
    rb.pose.position.setZero();
    rb.pose.rotation.setZero();
    return rb;
}

RigidBody
create_body(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
{
    return create_body(vertices, edges, Eigen::MatrixXi());
}

TEST_CASE("Rigid edge-vertex time of impact", "[ccd][rigid_toi][edge_vertex]")
{
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

    Pose<double> bodyA_pose_t0(
        Eigen::Vector3d(0, 0.5, 0).head(pos_ndof),
        Eigen::Vector3d(0, 0, 0).head(rot_ndof));
    Pose<double> bodyB_pose_t0 = Pose<double>::Zero(dim);
    Pose<double> bodyA_pose_t1 = Pose<double>::Zero(dim);
    Pose<double> bodyB_pose_t1 = Pose<double>::Zero(dim);

    double expected_toi = -1;
    bool is_impact_expected = true;
    SECTION("Translation")
    {
        double posx = GENERATE(-2 - 1e-8, -2, -1, -0.5);
        double sign = GENERATE(-1, 1);
        bodyA_vertices.row(0).x() = sign * posx;
        bodyA_vertices.row(1).x() = -sign * posx;
        double y_t1 = GENERATE(-10, -0.5, -1e-12, 0.0, 1e-12, 0.5, 10);
        expected_toi =
            bodyA_pose_t0.position.y() / (-y_t1 + bodyA_pose_t0.position.y());
        bodyA_pose_t1.position = Eigen::Vector3d(0, y_t1, 0).head(pos_ndof);

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
        bodyA_pose_t1.position = bodyA_pose_t0.position;
        bodyA_pose_t1.rotation = Eigen::Vector3d(0, 0, theta).tail(rot_ndof);
    }

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    double toi;
    bool is_impacting = compute_edge_vertex_time_of_impact(
        bodyA, bodyA_pose_t0, bodyA_pose_t1, /*vertex_id=*/0, //
        bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edge_id=*/0,   //
        toi);
    CAPTURE(toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    // if (is_impacting) {
    //     CHECK(
    //         toi
    //         == Approx(expected_toi)
    //                .margin(Constants::INTERVAL_ROOT_FINDER_TOL));
    // }
}

TEST_CASE("Rigid edge-edge time of impact", "[ccd][rigid_toi][edge_edge]")
{
    int dim = 3;
    int ndof = Pose<double>::dim_to_ndof(dim);
    int pos_ndof = Pose<double>::dim_to_pos_ndof(dim);
    int rot_ndof = Pose<double>::dim_to_rot_ndof(dim);

    Eigen::MatrixXd bodyA_vertices(2, dim);
    bodyA_vertices.row(0) << -1, 0, 0;
    bodyA_vertices.row(1) << 1, 0, 0;
    Eigen::MatrixXd bodyB_vertices(2, dim);
    bodyB_vertices.row(0) << 0, 0, -1;
    bodyB_vertices.row(1) << 0, 0, 1;

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    bodyB_edges.row(0) << 0, 1;

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    Pose<double> bodyA_pose_t0 = Pose<double>::Zero(dim);
    double z = GENERATE(1000.0, 10.0, 1.0 + 1e-8, 1.0, 1.0 - 1e-8, 0.5, 0.0);
    double sign = GENERATE(-1, 1);
    double z_t0 = sign * z;
    bodyA_pose_t0.position.y() = 1.0;
    bodyA_pose_t0.position.z() = z_t0;
    Pose<double> bodyB_pose_t0 = Pose<double>::Zero(dim);
    Pose<double> bodyA_pose_t1 = Pose<double>::Zero(dim);
    bodyA_pose_t1.position.z() = z_t0;
    Pose<double> bodyB_pose_t1 = Pose<double>::Zero(dim);

    double expected_toi = -1;
    bool is_impact_expected = false;
    SECTION("Translation")
    {
        double y_t1 =
            GENERATE(2.0, 1.0, 0.1, 1e-6, 0.0, -1e-6, -1, -2, -10, -1000);
        expected_toi = 1 / (1 - y_t1);
        is_impact_expected = y_t1 <= 0.0 && z_t0 >= -1 && z_t0 <= 1;
        bodyA_pose_t1.position.y() = y_t1;
    }

    double toi;
    bool is_impacting = compute_edge_edge_time_of_impact(
        bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
        bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeA_id=*/0, //
        toi);
    CAPTURE(
        bodyA_pose_t0.position.transpose(), bodyA_pose_t1.position.transpose(),
        bodyA.world_vertex(bodyA_pose_t0, 0).transpose(),
        bodyA.world_vertex(bodyA_pose_t1, 0).transpose(), toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    if (is_impacting) {
        CHECK(
            toi
            == Approx(expected_toi).margin(Constants::SCREWING_CCD_LENGTH_TOL));
    }
}

TEST_CASE("Rigid face-vertex time of impact", "[ccd][rigid_toi][face_vertex]")
{
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

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges, bodyA_faces);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges, bodyB_faces);

    Pose<double> bodyA_pose_t0 = bodyA.pose;
    Pose<double> bodyB_pose_t0 = bodyB.pose;

    Pose<double> bodyA_pose_t1 = bodyA.pose;
    Pose<double> bodyB_pose_t1 = bodyB.pose;

    double expected_toi = -1;
    bool is_impact_expected = false;
    SECTION("Translation")
    {
        double z_t1 = GENERATE(1, 1e-8, 0.0, -1e-8, -2.0, -10);
        expected_toi = 1 / (-z_t1 + 1);
        is_impact_expected = z_t1 <= 0.0 && y >= 0 && y <= 1.0;
        bodyB_pose_t1.position.z() = z_t1 - bodyB.vertices(0, 2);
    }
    // SECTION("Rotation")
    // {
    //     expected_toi = ;
    //     is_impact_expected = ;
    //     bodyB_velocity.rotation.z() = ;
    // }

    CAPTURE(
        y, bodyB.vertices.row(0), bodyB_pose_t0.position.transpose(),
        bodyB_pose_t1.position.transpose(),
        bodyB.world_vertex(bodyB_pose_t0, 0).transpose(),
        bodyB.world_vertex(bodyB_pose_t1, 0).transpose());
    double toi;
    bool is_impacting = compute_face_vertex_time_of_impact(
        bodyB, bodyB_pose_t0, bodyB_pose_t1, /*vertex_id=*/0, // Vertex body
        bodyA, bodyA_pose_t0, bodyA_pose_t1, /*face_id=*/0,   // Face body
        // Output time of impact
        toi);
    CAPTURE(toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    if (is_impacting) {
        CHECK(
            toi
            == Approx(expected_toi).margin(Constants::SCREWING_CCD_LENGTH_TOL));
    }
}

TEST_CASE("Fast EE case", "[!benchmark][ccd][rigid_toi][edge_edge][fast]")
{
    Eigen::MatrixXd bodyA_vertices = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd bodyB_vertices = Eigen::MatrixXd::Zero(2, 3);

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    bodyB_edges.row(0) << 0, 1;

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    // clang-format off
    bodyA.vertices.row(0) << 2.66473512640082, 0.622074238426736, 0.0506824409538513;
    bodyA.vertices.row(1) << 2.11663114018755, 0.24694623070118, 0.689392835886464;

    bodyB.vertices.row(0) << -8.23540431483827, -2.00204356583054, -0.0850398676470792;
    bodyB.vertices.row(1) << -7.59698634143846, -2.65778542381018, 0.220598272033643;

    Pose<double> bodyA_pose_t0(
        Eigen::Vector3d(-0.0743518262648221, 2.0045466941596, -7.30362621222097e-05),
        Eigen::Vector3d(1.56635729875648, 0.00484560212950007, -0.00477941835507233)
    );
    Pose<double> bodyA_pose_t1(
        Eigen::Vector3d(-0.0743518262648221, 2.0045466941596, -7.30362621222097e-05),
        Eigen::Vector3d(1.56635729875648, 0.00484560212950007, -0.00477941835507233)
    );
    Pose<double> bodyB_pose_t0(
        Eigen::Vector3d(-0.0625913211582116, 9.56420995832003, -0.0276191274028154),
        Eigen::Vector3d(-0.0330556101141157, -0.0327790458706025, 1.56913520538901)
    );
    Pose<double> bodyB_pose_t1(
        Eigen::Vector3d(-0.0625913211582116, 9.56371945832003, -0.0276191274028154),
        Eigen::Vector3d(-0.0467509876415133, -0.0464937186371367, 1.56888122389168)
    );
    // clang-format on

    BENCHMARK("Fast EE case")
    {
        double toi;
        bool is_impacting = compute_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeA_id=*/0, //
            toi);
    };
}

TEST_CASE("Slow EE case", "[!benchmark][ccd][rigid_toi][edge_edge][slow]")
{
    Eigen::MatrixXd bodyA_vertices = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd bodyB_vertices = Eigen::MatrixXd::Zero(2, 3);

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    bodyB_edges.row(0) << 0, 1;

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    // clang-format off
    bodyA.vertices.row(0) << -1.45054325721069, 2.29538642849017, 0.461193969193908;
    bodyA.vertices.row(1) << -1.12537372801783, 1.78188213954136, 1.12303701209507;

    bodyB.vertices.row(0) << -8.63773122773594, 0.523601125992661, 1.91725528075351;
    bodyB.vertices.row(1) << -8.00096703934998, 1.01814387645424, 2.44728456523133;

    Pose<double> bodyA_pose_t0(
        Eigen::Vector3d(-0.0743518262648221, 2.0045466941596, -7.30362621222097e-05),
        Eigen::Vector3d(1.56635729875648, 0.00484560212950006, -0.00477941835507233)
    );
    Pose<double> bodyA_pose_t1(
        Eigen::Vector3d(-0.0743518262648221, 2.0045466941596, -7.30362621222097e-05),
        Eigen::Vector3d(1.56635729875648, 0.00484560212950006, -0.00477941835507233)
    );
    Pose<double> bodyB_pose_t0(
        Eigen::Vector3d(-0.0625913211582116, 9.53281795832003, -0.0276191274028154),
        Eigen::Vector3d(-0.142589758121075, -0.142468016151593, 1.56466991422561)
    );
    Pose<double> bodyB_pose_t1(
        Eigen::Vector3d(-0.0625913211582116, 9.52447945832003, -0.0276191274028154),
        Eigen::Vector3d(-0.156274624583429, -0.156172306264317, 1.56372038452934)
    );
    // clang-format on

    BENCHMARK("Slow EE Case")
    {
        double toi;
        bool is_impacting = compute_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeA_id=*/0, //
            toi);
    };
}

TEST_CASE("Actual EE Collision", "[!benchmark][ccd][rigid_toi][edge_edge]")
{
    Eigen::MatrixXd bodyA_vertices = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd bodyB_vertices = Eigen::MatrixXd::Zero(2, 3);

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    bodyB_edges.row(0) << 0, 1;

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    // clang-format off
    bodyA.vertices.row(0) << -1, 0, 0;
    bodyA.vertices.row(1) << 1, 0, 0;

    bodyB.vertices.row(0) << 0, 0, -1;
    bodyB.vertices.row(1) << 0, 0, 1;

    Pose<double> bodyA_pose_t0(
        Eigen::Vector3d(0, 0.5, 0),
        Eigen::Vector3d(0, 0, 0)
    );
    Pose<double> bodyA_pose_t1(
        Eigen::Vector3d(0, 0.5, 0),
        Eigen::Vector3d(0, 0, 0)
    );
    Pose<double> bodyB_pose_t0(
        Eigen::Vector3d(0, 1, 0),
        Eigen::Vector3d(0, 0, 0)
    );
    Pose<double> bodyB_pose_t1(
        Eigen::Vector3d(0, 1, 0),
        Eigen::Vector3d(91 * 3.14, 0, 0)
    );
    // clang-format on

    BENCHMARK("Actually EE Collision")
    {
        double toi;
        bool is_impacting = compute_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeA_id=*/0, //
            toi);
    };
}
