#include <catch2/catch.hpp>

#include <boost/filesystem>
#include <string>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>
#include <nlohmann/json.hpp>

using namespace ccd;
using namespace nlohmann;

TEST_CASE("AABB initilization", "[hashgrid][AABB]")
{
    int dim = GENERATE(2, 3);
    CAPTURE(dim);
    AABB aabb;
    Eigen::VectorXd actual_center(dim);
    SECTION("Empty AABB")
    {
        aabb = AABB(Eigen::VectorXd::Zero(dim), Eigen::VectorXd::Zero(dim));
        actual_center = Eigen::VectorXd::Zero(dim);
    }
    SECTION("Box centered at zero")
    {
        Eigen::VectorXd min =
            Eigen::VectorXd::Random(dim).array() - 1; // in range [-2, 0]
        Eigen::VectorXd max = -min;
        aabb = AABB(min, max);
        actual_center = Eigen::VectorXd::Zero(dim);
    }
    SECTION("Box not centered at zero")
    {
        Eigen::VectorXd min(dim), max(dim);
        if (dim == 2) {
            min << 5.1, 3.14;
            max << 10.4, 7.89;
            actual_center << 7.75, 5.515;
        } else {
            min << 5.1, 3.14, 7.94;
            max << 10.4, 7.89, 10.89;
            actual_center << 7.75, 5.515, 9.415;
        }
        aabb = AABB(min, max);
    }
    Eigen::VectorXd center_diff = aabb.getCenter() - actual_center;
    CHECK(center_diff.norm() == Approx(0.0).margin(1e-12));
}

TEST_CASE("AABB overlapping", "[hashgrid][AABB]")
{
    AABB a, b;
    bool are_overlaping = false;
    SECTION("a to the right of b")
    {
        a = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlaping = false;
        }
    }
    SECTION("b to the right of a")
    {
        b = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlaping = false;
        }
    }
    SECTION("a above b")
    {
        a = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlaping = false;
        }
    }
    SECTION("a above b")
    {
        b = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlaping = false;
        }
    }
    CHECK(AABB::are_overlaping(a, b) == are_overlaping);
}

TEST_CASE("2D hash grid", "[hashgrid][2D]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd displacements;
    SECTION("Simple")
    {
        vertices.resize(10, 2);
        vertices << -2.4375, 0.0, 0.0, 0.0, 2.4375, 0.0, -3.0, -3.0, 0.0, -3.0,
            3.0, -3.0, -3.0, -6.0, -0.5625, -6.0, 0.5625, -6.0, 3.0, -6.0;

        edges.resize(9, 2);
        edges << 0, 1, 1, 2, 1, 4, 3, 4, 4, 5, 3, 6, 5, 9, 6, 7, 8, 9;
    }
    SECTION("Complex")
    {
        vertices.resize(64, 2);
        vertices << -10.0, -1.0, 10.0, -1.0, 10.0, 0.0, -10.0, 0.0,
            -12.260907953111856, -5.401392778833883, -13.256614411075963,
            -5.308825778404078, -13.34918141150577, -6.304532236368185,
            -12.353474953541662, -6.39709923679799, -10.205723196569847,
            -1.650089676888279, -11.19710985826582, -1.5191221732420903,
            -11.32807736191201, -2.5105088349380633, -10.336690700216037,
            -2.641476338584252, -10.233710759216248, 1.0244045371628476,
            -10.252058690562622, 0.024572874623943697, -9.252227028023718,
            0.00622494327756784, -9.233879096677343, 1.0060566058164717,
            12.30976068250293, -7.0522114282283965, 11.370487560538956,
            -7.395381943214828, 11.713658075525387, -8.3346550651788,
            12.65293119748936, -7.99148455019237, 15.358285786490363,
            -12.974846537789217, 14.367178024999136, -13.107908196849401,
            14.50023968405932, -14.099015958340628, 15.491347445550547,
            -13.965954299280444, -12.287270890151099, 0.6284868305968881,
            -11.846578218373745, -0.26917121623483126, -10.948920171542024,
            0.17152145554252315, -11.389612843319378, 1.0691795023742425,
            -4.160227868182811, 0.00017661538729152326, -3.160877050978387,
            0.036203607391201265, -3.1969040429822964, 1.0355544245956256,
            -4.196254860186721, 0.9995274325917158, 10.153702665044973,
            -0.002568366598117877, 10.17203116561408, 0.9972636523265054,
            9.172199146689456, 1.0155921528956124, 9.15387064612035,
            0.015760133970989076, 10.632702731376751, -2.461378460585955,
            11.625885082675442, -2.5779495422005905, 11.742456164290077,
            -1.5847671909019003, 10.749273812991387, -1.468196109287265,
            -9.189920836698171, 0.0046536337329208255, -8.18992896640201,
            0.0006213463443575096, -8.185896679013446, 1.0006132166405195,
            -9.185888549309608, 1.004645504029083, 1.3030037205605254,
            0.01235562685364272, 2.302973744436933, 0.004612801231432595,
            2.310716570059143, 1.0045828251078406, 1.3107465461827355,
            1.0123256507300507, 9.149900842491041, 0.007083139789913373,
            9.164573010886311, 1.0069754977337977, 8.164680652942428,
            1.0216476661290685, 8.150008484547158, 0.021755308185184286,
            -5.7155892517981455, 1.21614306103357e-05, -4.718344360077054,
            0.07419184824706188, -4.792524046893505, 1.0714367399681546,
            -5.789768938614597, 0.9972570531517031, -0.5052156587592828,
            0.024789052181432847, 0.4945801789736279, 0.004583041888764683,
            0.514786189266296, 1.0043788796216757, -0.4850096484666146,
            1.024584889914344, -6.955795314526457, 0.005049336769379409,
            -6.960731719921333, 1.0050371526460418, -7.960719535797995,
            1.0001007472511645, -7.955783130403119, 0.00011293137450218982;

        edges.resize(64, 2);
        edges << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 8, 9, 9, 10,
            10, 11, 11, 8, 12, 13, 13, 14, 14, 15, 15, 12, 16, 17, 17, 18, 18,
            19, 19, 16, 20, 21, 21, 22, 22, 23, 23, 20, 24, 25, 25, 26, 26, 27,
            27, 24, 28, 29, 29, 30, 30, 31, 31, 28, 32, 33, 33, 34, 34, 35, 35,
            32, 36, 37, 37, 38, 38, 39, 39, 36, 40, 41, 41, 42, 42, 43, 43, 40,
            44, 45, 45, 46, 46, 47, 47, 44, 48, 49, 49, 50, 50, 51, 51, 48, 52,
            53, 53, 54, 54, 55, 55, 52, 56, 57, 57, 58, 58, 59, 59, 56, 60, 61,
            61, 62, 62, 63, 63, 60;
    }

    displacements.resize(vertices.rows(), vertices.cols());
    for (int i = 0; i < 10; i++) {
        displacements.setRandom();
        displacements *= 10;

        EdgeVertexImpacts brute_force_ev_impacts;
        detect_edge_vertex_collisions(
            vertices, displacements, edges, brute_force_ev_impacts,
            DetectionMethod::BRUTE_FORCE);

        EdgeVertexImpacts hash_ev_impacts;
        detect_edge_vertex_collisions(
            vertices, displacements, edges, hash_ev_impacts,
            DetectionMethod::HASH_GRID);

        REQUIRE(brute_force_ev_impacts.size() == hash_ev_impacts.size());
        std::sort(
            brute_force_ev_impacts.begin(), brute_force_ev_impacts.end(),
            compare_impacts_by_time<EdgeVertexImpact>);
        std::sort(
            hash_ev_impacts.begin(), hash_ev_impacts.end(),
            compare_impacts_by_time<EdgeVertexImpact>);
        CHECK(brute_force_ev_impacts == hash_ev_impacts);
    }
}

TEST_CASE("3D hash grid", "[hashgrid][3D]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXd displacements;
    Eigen::MatrixXi edges;
    Eigen::MatrixXi faces;
    Eigen::VectorXi group_ids;

    SECTION("Simple")
    {
        vertices.resize(4, 3);
        vertices.row(0) << -1, -1, 0;
        vertices.row(1) << 1, -1, 0;
        vertices.row(2) << 0, 1, 1;
        vertices.row(3) << 0, 1, -1;

        edges.resize(2, 2);
        edges.row(0) << 0, 1;
        edges.row(1) << 2, 3;

        SECTION("Without group ids") {}
        SECTION("With group ids")
        {
            group_ids.resize(4);
            group_ids << 0, 0, 1, 1;
        }

        faces.resize(0, 3);

        displacements = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        displacements.col(1).head(2).setConstant(2);
        displacements.col(1).tail(2).setConstant(-2);
    }
    SECTION("Complex")
    {
        std::string fname = GENERATE(std::string("cube.obj")
                                     /*, std::string("bunny-lowpoly.obj")*/);

        boost::filesystem::path mesh_path = std::filesystem::path(__FILE__)
                                                .parent_path()
                                                .parent_path()
                                                .parent_path()
            / "meshes" / fname;
        igl::read_triangle_mesh(mesh_path.string(), vertices, faces);
        igl::edges(faces, edges);

        displacements = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        displacements.col(1).setOnes();
    }

    for (int i = 0; i < 10; i++) {
        EdgeVertexImpacts brute_force_ev_impacts;
        EdgeEdgeImpacts brute_force_ee_impacts;
        FaceVertexImpacts brute_force_fv_impacts;
        detect_collisions(
            vertices, displacements, edges, faces, group_ids,
            CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
            brute_force_ev_impacts, brute_force_ee_impacts,
            brute_force_fv_impacts, DetectionMethod::BRUTE_FORCE);
        REQUIRE(brute_force_ev_impacts.size() == 0);

        EdgeVertexImpacts hash_ev_impacts;
        EdgeEdgeImpacts hash_ee_impacts;
        FaceVertexImpacts hash_fv_impacts;
        detect_collisions(
            vertices, displacements, edges, faces, group_ids,
            CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
            hash_ev_impacts, hash_ee_impacts, hash_fv_impacts,
            DetectionMethod::HASH_GRID);
        REQUIRE(hash_ev_impacts.size() == 0);

        REQUIRE(brute_force_ee_impacts.size() == hash_ee_impacts.size());
        std::sort(
            brute_force_ee_impacts.begin(), brute_force_ee_impacts.end(),
            compare_impacts_by_time<EdgeEdgeImpact>);
        std::sort(
            hash_ee_impacts.begin(), hash_ee_impacts.end(),
            compare_impacts_by_time<EdgeEdgeImpact>);
        CHECK(brute_force_ee_impacts == hash_ee_impacts);

        REQUIRE(brute_force_fv_impacts.size() == hash_fv_impacts.size());
        std::sort(
            brute_force_fv_impacts.begin(), brute_force_fv_impacts.end(),
            compare_impacts_by_time<FaceVertexImpact>);
        std::sort(
            hash_fv_impacts.begin(), hash_fv_impacts.end(),
            compare_impacts_by_time<FaceVertexImpact>);
        CHECK(brute_force_fv_impacts == hash_fv_impacts);

        displacements.setRandom();
        displacements *= 3;
    }
}
