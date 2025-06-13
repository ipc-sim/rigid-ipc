#include <catch2/catch.hpp>

#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <ccd/ccd.hpp>
#include <logger.hpp>

using namespace ipc::rigid;

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

        Impacts brute_force_impacts;
        detect_collisions(
            vertices, vertices + displacements, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), CollisionType::EDGE_VERTEX, brute_force_impacts,
            DetectionMethod::BRUTE_FORCE);

        Impacts hash_impacts;
        detect_collisions(
            vertices, vertices + displacements, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), CollisionType::EDGE_VERTEX, hash_impacts,
            DetectionMethod::HASH_GRID);

        REQUIRE(brute_force_impacts.size() == hash_impacts.size());
        std::sort(
            brute_force_impacts.ev_impacts.begin(),
            brute_force_impacts.ev_impacts.end(),
            compare_impacts_by_time<EdgeVertexImpact>);
        std::sort(
            hash_impacts.ev_impacts.begin(), hash_impacts.ev_impacts.end(),
            compare_impacts_by_time<EdgeVertexImpact>);
        CHECK(brute_force_impacts.ev_impacts == hash_impacts.ev_impacts);
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

        SECTION("Without group ids") { }
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
        std::string fname = GENERATE(
            std::string("cube.obj")
            /*, std::string("bunny-lowpoly.obj")*/);

        fs::path mesh_path =
            fs::path(__FILE__).parent_path().parent_path().parent_path()
            / "meshes" / fname;
        igl::read_triangle_mesh(mesh_path.string(), vertices, faces);
        igl::edges(faces, edges);

        displacements = Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
        displacements.col(1).setOnes();
    }

    for (int i = 0; i < 10; i++) {
        Impacts brute_force_impacts;
        detect_collisions(
            vertices, vertices + displacements, edges, faces, group_ids,
            CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
            brute_force_impacts, DetectionMethod::BRUTE_FORCE);
        REQUIRE(brute_force_impacts.ev_impacts.size() == 0);

        Impacts hash_impacts;
        detect_collisions(
            vertices, vertices + displacements, edges, faces, group_ids,
            CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX, hash_impacts,
            DetectionMethod::HASH_GRID);
        REQUIRE(hash_impacts.ev_impacts.size() == 0);

        CAPTURE(i);
        if (brute_force_impacts.ee_impacts.size()
            != hash_impacts.ee_impacts.size()) {
            std::cout << displacements << std::endl;
        }
        REQUIRE(
            brute_force_impacts.ee_impacts.size()
            == hash_impacts.ee_impacts.size());
        std::sort(
            brute_force_impacts.ee_impacts.begin(),
            brute_force_impacts.ee_impacts.end(),
            compare_impacts_by_time<EdgeEdgeImpact>);
        std::sort(
            hash_impacts.ee_impacts.begin(), hash_impacts.ee_impacts.end(),
            compare_impacts_by_time<EdgeEdgeImpact>);
        bool is_equal =
            brute_force_impacts.ee_impacts == hash_impacts.ee_impacts;
        if (!is_equal && brute_force_impacts.size() > 0) {
            spdlog::error(
                "bf_impacts.ee[0]={{time={:g}, e0i={:d}, α₀={:g}, e0i={:d}, "
                "α₁={:g}}}",
                brute_force_impacts.ee_impacts[0].time,
                brute_force_impacts.ee_impacts[0].impacted_edge_index,
                brute_force_impacts.ee_impacts[0].impacted_alpha,
                brute_force_impacts.ee_impacts[0].impacting_edge_index,
                brute_force_impacts.ee_impacts[0].impacting_alpha);
            spdlog::error(
                "hash_impacts[0]={{time={:g}, e0i={:d}, α₀={:g}, e0i={:d}, "
                "α₁={:g}}}",
                hash_impacts.ee_impacts[0].time,
                hash_impacts.ee_impacts[0].impacted_edge_index,
                hash_impacts.ee_impacts[0].impacted_alpha,
                hash_impacts.ee_impacts[0].impacting_edge_index,
                hash_impacts.ee_impacts[0].impacting_alpha);
        }
        CHECK(is_equal);

        if (brute_force_impacts.fv_impacts.size()
            != hash_impacts.fv_impacts.size()) {
            std::cout << logger::fmt_eigen(vertices) << std::endl;
            std::cout << logger::fmt_eigen(displacements) << std::endl;
            std::cout << edges << std::endl;
            std::cout << faces << std::endl;
            std::cout << group_ids << std::endl;
        }
        REQUIRE(
            brute_force_impacts.fv_impacts.size()
            == hash_impacts.fv_impacts.size());
        std::sort(
            brute_force_impacts.fv_impacts.begin(),
            brute_force_impacts.fv_impacts.end(),
            compare_impacts_by_time<FaceVertexImpact>);
        std::sort(
            hash_impacts.fv_impacts.begin(), hash_impacts.fv_impacts.end(),
            compare_impacts_by_time<FaceVertexImpact>);
        CHECK(brute_force_impacts.fv_impacts == hash_impacts.fv_impacts);

        displacements.setRandom();
        displacements *= 3;
    }
}

TEST_CASE("3D hash grid case 1", "[hashgrid][3D]")
{
    // clang-format off
    Eigen::MatrixXd vertices(8, 3);
    vertices <<
        -0.5, -0.5, -0.5,
         0.5,  0.5, -0.5,
         0.5, -0.5, -0.5,
        -0.5,  0.5, -0.5,
        -0.5,  0.5,  0.5,
        -0.5, -0.5,  0.5,
         0.5,  0.5,  0.5,
         0.5, -0.5,  0.5;
    Eigen::MatrixXd displacements(8, 3);
    displacements <<
        -0.521693493948175,  2.8007648996081,   -1.03474376072862,
        -2.10255278698288,  2.45566771340355,  -2.93838656597696,
         2.39530917880838, -1.59274082658474,   0.53698562529729,
        -2.03863176751818,  2.80492759021229,   1.11740437155469,
         2.71588332192781,  0.418008697879505,  0.215272719606419,
        -2.14900835936377, -0.527814739163879,  0.0885984250756906,
         1.6165041730816, -2.98232112730961,   1.07373024712956,
         0.585636982501315,  0.128813307326666, -1.81573649347561;
    Eigen::MatrixXi edges(18, 2);
    edges <<
        0, 1,
        0, 2,
        1, 2,
        0, 3,
        1, 3,
        0, 4,
        3, 4,
        0, 5,
        4, 5,
        1, 6,
        2, 6,
        3, 6,
        4, 6,
        5, 6,
        0, 7,
        2, 7,
        5, 7,
        6, 7;
    Eigen::MatrixXi faces(12, 3);
    faces <<
        0, 1, 2,
        0, 3, 1,
        0, 4, 3,
        0, 5, 4,
        3, 6, 1,
        3, 4, 6,
        2, 1, 6,
        2, 6, 7,
        0, 2, 7,
        0, 7, 5,
        5, 7, 6,
        5, 6, 4;
    // clang-format on
    Eigen::VectorXi group_ids;

    Impacts bf_impacts;
    detect_collisions(
        vertices, vertices + displacements, edges, faces, group_ids,
        CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX, bf_impacts,
        DetectionMethod::BRUTE_FORCE);
    REQUIRE(bf_impacts.ev_impacts.size() == 0);

    Impacts hash_impacts;
    detect_collisions(
        vertices, vertices + displacements, edges, faces, group_ids,
        CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX, hash_impacts,
        DetectionMethod::HASH_GRID);
    REQUIRE(hash_impacts.ev_impacts.size() == 0);

    REQUIRE(bf_impacts.ee_impacts.size() == hash_impacts.ee_impacts.size());
    std::sort(
        bf_impacts.ee_impacts.begin(), bf_impacts.ee_impacts.end(),
        compare_impacts_by_time<EdgeEdgeImpact>);
    std::sort(
        hash_impacts.ee_impacts.begin(), hash_impacts.ee_impacts.end(),
        compare_impacts_by_time<EdgeEdgeImpact>);
    bool is_equal = bf_impacts.ee_impacts == hash_impacts.ee_impacts;
    CHECK(is_equal);

    // for (int i = 0; i < bf_impacts.fv_impacts.size(); i++) {
    //     fmt::print(
    //         "bf_impacts.fv[{:d}]={{toi={:g}, fi={:d}, u={:g}, v={:g}, "
    //         "vi={:d}}}\n",
    //         i, bf_impacts.fv_impacts[i].time,
    //         bf_impacts.fv_impacts[i].face_index, bf_impacts.fv_impacts[i].u,
    //         bf_impacts.fv_impacts[i].v,
    //         bf_impacts.fv_impacts[i].vertex_index);
    // }
    // for (int i = 0; i < hash_impacts.fv_impacts.size(); i++) {
    //     fmt::print(
    //         "hash_impacts.fv[{:d}]={{toi={:g}, fi={:d}, u={:g}, v={:g}, "
    //         "vi={:d}}}\n",
    //         i, hash_impacts.fv_impacts[i].time,
    //         hash_impacts.fv_impacts[i].face_index, //
    //         hash_impacts.fv_impacts[i].u, hash_impacts.fv_impacts[i].v,
    //         hash_impacts.fv_impacts[i].vertex_index);
    // }
    // ipc::Candidates hash_candidates;
    // detect_collision_candidates_hash_grid(
    //     vertices, vertices + displacements, edges, faces, group_ids,
    //     CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
    //     hash_candidates);
    // for (int i = 0; i < hash_candidates.fv_candidates.size(); i++) {
    //     fmt::print(
    //         "hash_candidates.fv[{:d}]={{fi={:d}, vi={:d}}}\n", i,
    //         hash_candidates.fv_candidates[i].face_index,
    //         hash_candidates.fv_candidates[i].vertex_index);
    // }

    REQUIRE(bf_impacts.fv_impacts.size() == hash_impacts.fv_impacts.size());
    std::sort(
        bf_impacts.fv_impacts.begin(), bf_impacts.fv_impacts.end(),
        compare_impacts_by_time<FaceVertexImpact>);
    std::sort(
        hash_impacts.fv_impacts.begin(), hash_impacts.fv_impacts.end(),
        compare_impacts_by_time<FaceVertexImpact>);
    CHECK(bf_impacts.fv_impacts == hash_impacts.fv_impacts);
}

TEST_CASE("3D brute force is duplicate free", "[ccd][brute_force]")
{
    Eigen::MatrixXd vertices;
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

        SECTION("Without group ids") { }
        SECTION("With group ids")
        {
            group_ids.resize(4);
            group_ids << 0, 0, 1, 1;
        }

        faces.resize(0, 3);
    }
    SECTION("Complex")
    {
        std::string fname = GENERATE(
            std::string("cube.obj")
            /*, std::string("bunny-lowpoly.obj")*/);

        fs::path mesh_path =
            fs::path(__FILE__).parent_path().parent_path().parent_path()
            / "meshes" / fname;
        igl::read_triangle_mesh(mesh_path.string(), vertices, faces);
        igl::edges(faces, edges);
    }

    for (int i = 0; i < 10; i++) {
        ipc::Candidates candidates;
        using namespace CollisionType;
        detect_collision_candidates_brute_force(
            vertices, edges, faces, group_ids,
            (EDGE_VERTEX | EDGE_EDGE | FACE_VERTEX), candidates);

        tbb::parallel_sort(
            candidates.ev_candidates.begin(), candidates.ev_candidates.end());
        auto ev_unique_end = std::unique(
            candidates.ev_candidates.begin(), candidates.ev_candidates.end());
        CHECK(ev_unique_end == candidates.ev_candidates.end());
        CHECK(
            candidates.ev_candidates.size() <= edges.rows() * vertices.rows());

        tbb::parallel_sort(
            candidates.ee_candidates.begin(), candidates.ee_candidates.end());
        auto ee_unique_end = std::unique(
            candidates.ee_candidates.begin(), candidates.ee_candidates.end());
        CHECK(ee_unique_end == candidates.ee_candidates.end());
        CHECK(
            candidates.fv_candidates.size()
            <= edges.rows() * (edges.rows() - 1) / 2);

        tbb::parallel_sort(
            candidates.fv_candidates.begin(), candidates.fv_candidates.end());
        auto fv_unique_end = std::unique(
            candidates.fv_candidates.begin(), candidates.fv_candidates.end());
        CHECK(fv_unique_end == candidates.fv_candidates.end());
        CHECK(
            candidates.fv_candidates.size() <= faces.rows() * vertices.rows());
    }
}
