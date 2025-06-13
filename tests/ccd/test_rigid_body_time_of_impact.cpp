#include <catch2/catch.hpp>

#include <filesystem>
namespace fs = std::filesystem;

#include <igl/PI.h>
#include <igl/edges.h>

#include <ipc/distance/edge_edge.hpp>

// #include <ccd.hpp>
#include <ccd/piecewise_linear/time_of_impact.hpp>
#include <ccd/rigid/time_of_impact.hpp>
#include <constants.hpp>
#include <io/serialize_json.hpp>

using namespace ipc;
using namespace ipc::rigid;

const double TESTING_TOI_TOLERANCE = 1e-6;

void print_EE_obj(
    const RigidBody& bodyA,
    const Pose<double>& bodyA_pose_t0,
    const Pose<double>& bodyA_pose_t1,
    int edgeA_id,
    const RigidBody& bodyB,
    const Pose<double>& bodyB_pose_t0,
    const Pose<double>& bodyB_pose_t1,
    int edgeB_id,
    int n = 100)
{
    fmt::print("# Edge 1 vertices\n");
    for (int i = 0; i < n + 1; i++) {
        Pose<double> pose =
            Pose<double>::interpolate(bodyA_pose_t0, bodyA_pose_t1, i / n);
        std::cout
            << "v "
            << bodyA.world_vertex(pose, bodyA.edges(edgeA_id, 0)).transpose()
            << std::endl;
        std::cout
            << "v "
            << bodyA.world_vertex(pose, bodyA.edges(edgeA_id, 1)).transpose()
            << std::endl;
    }
    fmt::print("# Edge 2 vertices\n");
    for (int i = 0; i < n + 1; i++) {
        Pose<double> pose =
            Pose<double>::interpolate(bodyB_pose_t0, bodyB_pose_t1, i / n);
        std::cout
            << "v "
            << bodyB.world_vertex(pose, bodyB.edges(edgeB_id, 0)).transpose()
            << std::endl;
        std::cout
            << "v "
            << bodyB.world_vertex(pose, bodyB.edges(edgeB_id, 1)).transpose()
            << std::endl;
    }
    fmt::print("# Edge 1 surface\n");
    for (int i = 0; i < 4 * n + 2; i += 2) {
        if (i == 2 * n) {
            fmt::print("# Edge 2 surface\n");
            continue;
        }
        fmt::print("f {} {} {} {}\n", i + 1, i + 2, i + 4, i + 3);
    }
}

RigidBody create_body(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
{
    static int id = 0;
    int dim = vertices.cols();
    Pose<double> pose = Pose<double>::Zero(dim);
    RigidBody rb = RigidBody(
        vertices, edges, faces, pose,
        /*velocity=*/Pose<double>::Zero(pose.dim()),
        /*force=*/Pose<double>::Zero(pose.dim()),
        /*denisty=*/1.0,
        /*is_dof_fixed=*/VectorMax6b::Zero(pose.ndof()),
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
        toi, /*earliest_toi=*/1, /*toi_tolerance=*/TESTING_TOI_TOLERANCE);
    CAPTURE(toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    if (is_impacting) {
        // clang-format off
        CHECK(toi == Approx(expected_toi).margin(TESTING_TOI_TOLERANCE));
        // clang-format on
        CHECK(toi <= expected_toi);
    }
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
        toi, /*earliest_toi=*/1, /*toi_tolerance=*/TESTING_TOI_TOLERANCE);
    CAPTURE(
        bodyA_pose_t0.position.transpose(), bodyA_pose_t1.position.transpose(),
        bodyA.world_vertex(bodyA_pose_t0, 0).transpose(),
        bodyA.world_vertex(bodyA_pose_t1, 0).transpose(), toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    if (is_impacting) {
        // clang-format off
        CHECK(toi == Approx(expected_toi).margin(
            Constants::RIGID_CCD_LENGTH_TOL));
        // clang-format on
        CHECK(toi <= expected_toi);
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
        toi, /*earliest_toi=*/1, /*toi_tolerance=*/TESTING_TOI_TOLERANCE);
    CAPTURE(toi, expected_toi);
    CHECK(is_impacting == is_impact_expected);
    if (is_impact_expected) {
        // clang-format off
        CHECK(toi == Approx(expected_toi).margin(
            Constants::RIGID_CCD_LENGTH_TOL));
        // clang-format on
        CHECK(toi <= expected_toi);
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
        Eigen::Vector3d(2 * 3.14, 0, 0)
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

TEST_CASE("Actual VF Collision", "[!benchmark][ccd][rigid_toi][face_vertex]")
{
    Eigen::MatrixXd bodyA_vertices = Eigen::MatrixXd::Zero(3, 3);
    Eigen::MatrixXd bodyB_vertices = Eigen::MatrixXd::Zero(3, 3);

    Eigen::MatrixXi bodyA_faces(1, 3);
    bodyA_faces.row(0) << 0, 1, 2;
    Eigen::MatrixXi bodyB_faces(1, 3);
    bodyB_faces.row(0) << 0, 1, 2;

    Eigen::MatrixXi bodyA_edges, bodyB_edges;
    igl::edges(bodyA_faces, bodyA_edges);
    igl::edges(bodyB_faces, bodyB_edges);

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges, bodyA_faces);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges, bodyB_faces);

    // clang-format off
    bodyA.vertices.row(0) << 0.063161123153442, -0.00975209722618602, 0.0246948915619087;

    bodyB.vertices.row(0) << -0.00733894260009082, 0.0199670606490534, 0.000727755816038143;
    bodyB.vertices.row(1) << -0.0122514761614292, 0.0244249832266042, -0.00566776185395443;
    bodyB.vertices.row(2) << -0.00945035723923828, 0.025642170175986, -0.00591155842365654;

    Pose<double> bodyA_pose_t0(
        Eigen::Vector3d(-0.000591328227883731, 0.0888868556875028, -0.000277685809028307),
        Eigen::Vector3d(1.20966059163976, -1.2107299033619, -1.20913719904179)
    );
    Pose<double> bodyA_pose_t1(
        Eigen::Vector3d(-0.000591328227883731, 0.0879058556875029, -0.000277685809028307),
        Eigen::Vector3d(1.20966059163976, -1.2107299033619, -1.20913719904179)
    );
    Pose<double> bodyB_pose_t0(
        Eigen::Vector3d(-0.00052287436757153, 0.020051622753965, -4.06061775370095e-06),
        Eigen::Vector3d(1.21919966781284, -1.19787487380512, 1.19585545343065)
    );
    Pose<double> bodyB_pose_t1(
        Eigen::Vector3d(-0.00052287436757153, 0.020051622753965, -4.06061775370095e-06),
        Eigen::Vector3d(1.21919966781284, -1.19787487380512, 1.19585545343065)
    );
    // clang-format on

    BENCHMARK("Actually VF Collision")
    {
        double toi;
        bool is_impacting = compute_face_vertex_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*vertex_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*face_id=*/0,   //
            toi);
        // std::cout << toi << std::endl;
    };
}

TEST_CASE(
    "Extremly Slow EE Case",
    "[!benchmark][ccd][rigid_toi][edge_edge][extremly_slow]")
{
    Eigen::MatrixXd bodyA_vertices = Eigen::MatrixXd::Zero(2, 3);
    Eigen::MatrixXd bodyB_vertices = Eigen::MatrixXd::Zero(2, 3);

    Eigen::MatrixXi bodyA_edges(1, 2);
    bodyA_edges.row(0) << 0, 1;
    Eigen::MatrixXi bodyB_edges(1, 2);
    bodyB_edges.row(0) << 0, 1;

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges);

    Pose<double> bodyA_pose_t0, bodyA_pose_t1, bodyB_pose_t0, bodyB_pose_t1;

    double earliest_toi;

    SECTION("0")
    {
        // clang-format off
        bodyA.vertices.row(0) << 1.25, 0.625, -1.11022302462516e-16;
        bodyA.vertices.row(1) << 1.25, -0.625, -1.11022302462516e-16;

        bodyB.vertices.row(0) << 1.25, 0.625, 1.11022302462516e-16;
        bodyB.vertices.row(1) << 1.25, -0.625, 1.11022302462516e-16;

        bodyA_pose_t0 = Pose<double>(
            Eigen::Vector3d(-0.749789935368566, 1.00585262304029, 1.37760763963751e-05),
            Eigen::Vector3d(1.44479640128067, 0.727816581920491, 0.729072550558035)
        );
        bodyA_pose_t1 = Pose<double>(
            Eigen::Vector3d(-0.749767679498726, 0.999291290743525, 1.56515114692218e-05),
            Eigen::Vector3d(1.44714648553188, 0.720936564079057, 0.722362732765182)
        );
        bodyB_pose_t0 = Pose<double>(
            Eigen::Vector3d(0.749821260870553, 1.00573285487415, -1.25619191880717e-05),
            Eigen::Vector3d(0.836943953227619, 1.66096593034921, 1.66114035080701)
        );
        bodyB_pose_t1 = Pose<double>(
            Eigen::Vector3d(0.749804770233396, 0.999147283739376, -1.37411290439571e-05),
            Eigen::Vector3d(0.830693602998813, 1.6669129896604, 1.66713787700267)
        );
        // clang-format on
        earliest_toi = 0.739807;
    }
    SECTION("1")
    {
        // clang-format off
        bodyA.vertices.row(0) << 1.25, 0.625, -1.11022302462516e-16;
        bodyA.vertices.row(1) << 1.25, -0.625, -1.11022302462516e-16;

        bodyB.vertices.row(0) << 1.25, 0.625, 1.11022302462516e-16;
        bodyB.vertices.row(1) << 1.25, -0.625, 1.11022302462516e-16;

        bodyA_pose_t0 = Pose<double>(
            Eigen::Vector3d(-0.749781303981602, 1.00328869329824, 1.45053758187115e-05),
            Eigen::Vector3d(1.44571477658472, 0.725130939760901,
            0.726452486140648)
        );
        bodyA_pose_t1 = Pose<double>(
            Eigen::Vector3d(-0.749768505996024,
            0.999271155189446, 1.56109303515142e-05),
            Eigen::Vector3d(1.44714132170885, 0.720946434493765,
            0.722364003050481)
        );
        bodyB_pose_t0 = Pose<double>(
            Eigen::Vector3d(0.749814828894052, 1.0031128414029,
            -1.30219790400155e-05),
            Eigen::Vector3d(0.834519686991294, 1.66327347194605, 1.6634661858807)
        );
        bodyB_pose_t1 = Pose<double>(
            Eigen::Vector3d(0.749754900679728, 0.999118556601472,
            -1.68864775474044e-05),
            Eigen::Vector3d(0.830720747893929, 1.66690721127267, 1.66711165819808)
        );
        // clang-format on
        earliest_toi = 0.57421;
    }

    // print_EE_obj(
    //     bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
    //     bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeB_id=*/0);

    BENCHMARK("Extremly Slow EE CCD")
    {
        double toi;
        bool is_impacting = compute_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeB_id=*/0, //
            toi);
    };
}

TEST_CASE("Failing earliest tois", "[ccd][rigid_toi][failing_toi]")
{
    int id = GENERATE(range(0, 385));
    std::string filename =
        fmt::format("wrecking-ball/ccd-test-{:03d}.json", id);
    fs::path data_path =
        fs::path(__FILE__).parent_path().parent_path() / "data" / filename;

    std::ifstream input(data_path.string());
    nlohmann::json json = nlohmann::json::parse(input);

    Eigen::MatrixXd bodyA_vertices, bodyB_vertices;
    Eigen::MatrixXi bodyA_edges, bodyB_edges, bodyA_faces, bodyB_faces;
    Pose<double> bodyA_pose_t0, bodyA_pose_t1, bodyB_pose_t0, bodyB_pose_t1;

    std::string ccd_type = json["type"];

    nlohmann::json poseA_t0_json, poseA_t1_json, poseB_t0_json, poseB_t1_json;

    if (ccd_type == "ee") {
        Eigen::VectorXd tmp;

        bodyA_vertices.resize(2, 3);
        from_json(json["edge0"]["vertex0"], tmp);
        bodyA_vertices.row(0) = tmp;
        from_json(json["edge0"]["vertex1"], tmp);
        bodyA_vertices.row(1) = tmp;

        bodyB_vertices.resize(2, 3);
        from_json(json["edge1"]["vertex0"], tmp);
        bodyB_vertices.row(0) = tmp;
        from_json(json["edge1"]["vertex1"], tmp);
        bodyB_vertices.row(1) = tmp;

        bodyA_edges.resize(1, 2);
        bodyA_edges.row(0) << 0, 1;
        bodyB_edges.resize(1, 2);
        bodyB_edges.row(0) << 0, 1;

        poseA_t0_json = json["edge0"]["pose_t0"];
        poseA_t1_json = json["edge0"]["pose_t1"];
        poseB_t0_json = json["edge1"]["pose_t0"];
        poseB_t1_json = json["edge1"]["pose_t1"];
    } else if (ccd_type == "fv") {
        Eigen::VectorXd tmp;

        bodyA_vertices.resize(3, 3);
        from_json(json["face"]["vertex0"], tmp);
        bodyA_vertices.row(0) = tmp;
        from_json(json["face"]["vertex1"], tmp);
        bodyA_vertices.row(1) = tmp;
        from_json(json["face"]["vertex2"], tmp);
        bodyA_vertices.row(2) = tmp;

        bodyB_vertices.resize(1, 3);
        from_json(json["vertex"]["vertex"], tmp);
        bodyB_vertices.row(0) = tmp;

        bodyA_faces.resize(1, 3);
        bodyA_faces.row(0) << 0, 1, 2;
        bodyA_edges.resize(3, 2);
        bodyA_edges.row(0) << 0, 1;
        bodyA_edges.row(1) << 1, 2;
        bodyA_edges.row(2) << 2, 0;

        poseA_t0_json = json["face"]["pose_t0"];
        poseA_t1_json = json["face"]["pose_t1"];
        poseB_t0_json = json["vertex"]["pose_t0"];
        poseB_t1_json = json["vertex"]["pose_t1"];
    } else if (ccd_type == "ev") {
        // TODO
        return;
    }

    from_json(poseA_t0_json["position"], bodyA_pose_t0.position);
    from_json(poseA_t0_json["rotation"], bodyA_pose_t0.rotation);
    from_json(poseA_t1_json["position"], bodyA_pose_t1.position);
    from_json(poseA_t1_json["rotation"], bodyA_pose_t1.rotation);
    from_json(poseB_t0_json["position"], bodyB_pose_t0.position);
    from_json(poseB_t0_json["rotation"], bodyB_pose_t0.rotation);
    from_json(poseB_t1_json["position"], bodyB_pose_t1.position);
    from_json(poseB_t1_json["rotation"], bodyB_pose_t1.rotation);

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges, bodyA_faces);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges, bodyB_faces);

    double toi = 1;
    bool is_impacting = false;

    if (ccd_type == "ee") {
        is_impacting = compute_piecewise_linear_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeB_id=*/0, //
            toi);
    } else if (ccd_type == "fv") {
        is_impacting = compute_piecewise_linear_face_vertex_time_of_impact(
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*vertex_id=*/0, // Vertex body
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*face_id=*/0,   // Face body
            toi);
    }

    CAPTURE(ccd_type);
    if (!is_impacting) {
        return;
    }

    // toi *= 0.99;
    double toi2 = 1;

    Pose<double> bodyA_pose_toi =
        Pose<double>::interpolate(bodyA_pose_t0, bodyA_pose_t1, toi);
    Pose<double> bodyB_pose_toi =
        Pose<double>::interpolate(bodyB_pose_t0, bodyB_pose_t1, toi);

    if (ccd_type == "ee") {
        is_impacting = compute_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_toi, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_toi, /*edgeB_id=*/0, //
            toi2);
    } else if (ccd_type == "fv") {
        is_impacting = compute_face_vertex_time_of_impact(
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*vertex_id=*/0, // Vertex
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*face_id=*/0,   // Face
            toi2);
    }
    CAPTURE(toi, toi2, toi * toi2);
    CHECK(!is_impacting);
}

TEST_CASE("toi=0", "[ccd][rigid_toi][thisone]")
{
    int id = GENERATE(range(0, 13));
    std::string filename = fmt::format("kinematic/ccd-test-{:03d}.json", id);
    fs::path data_path =
        fs::path(__FILE__).parent_path().parent_path() / "data" / filename;

    std::ifstream input(data_path.string());
    nlohmann::json json = nlohmann::json::parse(input);

    Eigen::MatrixXd bodyA_vertices, bodyB_vertices;
    Eigen::MatrixXi bodyA_edges, bodyB_edges, bodyA_faces, bodyB_faces;
    Pose<double> bodyA_pose_t0, bodyA_pose_t1, bodyB_pose_t0, bodyB_pose_t1;

    std::string ccd_type = json["type"];

    nlohmann::json poseA_t0_json, poseA_t1_json, poseB_t0_json, poseB_t1_json;

    if (ccd_type == "ee") {
        Eigen::VectorXd tmp;

        bodyA_vertices.resize(2, 3);
        from_json(json["edge0"]["vertex0"], tmp);
        bodyA_vertices.row(0) = tmp;
        from_json(json["edge0"]["vertex1"], tmp);
        bodyA_vertices.row(1) = tmp;

        bodyB_vertices.resize(2, 3);
        from_json(json["edge1"]["vertex0"], tmp);
        bodyB_vertices.row(0) = tmp;
        from_json(json["edge1"]["vertex1"], tmp);
        bodyB_vertices.row(1) = tmp;

        bodyA_edges.resize(1, 2);
        bodyA_edges.row(0) << 0, 1;
        bodyB_edges.resize(1, 2);
        bodyB_edges.row(0) << 0, 1;

        poseA_t0_json = json["edge0"]["pose_t0"];
        poseA_t1_json = json["edge0"]["pose_t1"];
        poseB_t0_json = json["edge1"]["pose_t0"];
        poseB_t1_json = json["edge1"]["pose_t1"];
    } else if (ccd_type == "fv") {
        Eigen::VectorXd tmp;

        bodyA_vertices.resize(3, 3);
        from_json(json["face"]["vertex0"], tmp);
        bodyA_vertices.row(0) = tmp;
        from_json(json["face"]["vertex1"], tmp);
        bodyA_vertices.row(1) = tmp;
        from_json(json["face"]["vertex2"], tmp);
        bodyA_vertices.row(2) = tmp;

        bodyB_vertices.resize(1, 3);
        from_json(json["vertex"]["vertex"], tmp);
        bodyB_vertices.row(0) = tmp;

        bodyA_faces.resize(1, 3);
        bodyA_faces.row(0) << 0, 1, 2;
        bodyA_edges.resize(3, 2);
        bodyA_edges.row(0) << 0, 1;
        bodyA_edges.row(1) << 1, 2;
        bodyA_edges.row(2) << 2, 0;

        poseA_t0_json = json["face"]["pose_t0"];
        poseA_t1_json = json["face"]["pose_t1"];
        poseB_t0_json = json["vertex"]["pose_t0"];
        poseB_t1_json = json["vertex"]["pose_t1"];
    } else if (ccd_type == "ev") {
        // TODO
        return;
    }

    from_json(poseA_t0_json["position"], bodyA_pose_t0.position);
    from_json(poseA_t0_json["rotation"], bodyA_pose_t0.rotation);
    from_json(poseA_t1_json["position"], bodyA_pose_t1.position);
    from_json(poseA_t1_json["rotation"], bodyA_pose_t1.rotation);
    from_json(poseB_t0_json["position"], bodyB_pose_t0.position);
    from_json(poseB_t0_json["rotation"], bodyB_pose_t0.rotation);
    from_json(poseB_t1_json["position"], bodyB_pose_t1.position);
    from_json(poseB_t1_json["rotation"], bodyB_pose_t1.rotation);

    RigidBody bodyA = create_body(bodyA_vertices, bodyA_edges, bodyA_faces);
    RigidBody bodyB = create_body(bodyB_vertices, bodyB_edges, bodyB_faces);

    double toi = 1;
    bool is_impacting = false;

    if (ccd_type == "ee") {
        is_impacting = compute_piecewise_linear_edge_edge_time_of_impact(
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*edgeA_id=*/0, //
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*edgeB_id=*/0, //
            toi);
    } else if (ccd_type == "fv") {
        is_impacting = compute_piecewise_linear_face_vertex_time_of_impact(
            bodyB, bodyB_pose_t0, bodyB_pose_t1, /*vertex_id=*/0, // Vertex body
            bodyA, bodyA_pose_t0, bodyA_pose_t1, /*face_id=*/0,   // Face body
            toi);
    }

    Eigen::MatrixXd VA_t0 = bodyA.world_vertices(bodyA_pose_t0);
    Eigen::MatrixXd VB_t0 = bodyB.world_vertices(bodyB_pose_t0);

    if (ccd_type == "ee") {
        double distance = ipc::edge_edge_distance(
            VA_t0.row(0), VA_t0.row(1), VB_t0.row(0), VB_t0.row(1));
        // std::cout << "distance_t0=" << distance << std::endl;
    }

    Pose<double> bodyA_pose_toi =
        Pose<double>::interpolate(bodyA_pose_t0, bodyA_pose_t1, toi);
    Pose<double> bodyB_pose_toi =
        Pose<double>::interpolate(bodyB_pose_t0, bodyB_pose_t1, toi);
    Eigen::MatrixXd VA_toi = bodyA.world_vertices(bodyA_pose_toi);
    Eigen::MatrixXd VB_toi = bodyB.world_vertices(bodyB_pose_toi);

    if (ccd_type == "ee") {
        double distance = ipc::edge_edge_distance(
            VA_toi.row(0), VA_toi.row(1), VB_toi.row(0), VB_toi.row(1));
        // std::cout << "distance_toi=" << distance << std::endl;
        CHECK(distance == Approx(0).margin(1e-10));
    }

    CAPTURE(ccd_type);
    if (is_impacting) {
        CHECK(toi > 0);
    }
}
