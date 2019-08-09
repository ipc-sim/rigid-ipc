#include <catch2/catch.hpp>

#include <fstream>
#include <iostream>

#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>
#include <nlohmann/json.hpp>

using namespace ccd;
using namespace nlohmann;

TEST_CASE("AABB initilization", "[hashgrid][AABB]")
{
    AABB aabb;
    Eigen::Vector2d actual_center;
    SECTION("Empty AABB") { actual_center = Eigen::Vector2d::Zero(); }
    SECTION("Box centered at zero")
    {
        Eigen::Vector2d min
            = Eigen::Vector2d::Random().array() - 1; // in range [-2, 0]
        Eigen::Vector2d max = -min;
        aabb = AABB(min, max);
        actual_center = Eigen::Vector2d::Zero();
    }
    SECTION("Box not centered at zero")
    {
        Eigen::Vector2d min(5.1, 3.14), max(10.4, 7.89);
        aabb = AABB(min, max);
        actual_center = Eigen::Vector2d(7.75, 5.515);
    }
    Eigen::Vector2d center_diff = aabb.getCenter() - actual_center;
    CHECK(center_diff.x() == Approx(0.0).margin(1e-12));
    CHECK(center_diff.y() == Approx(0.0).margin(1e-12));
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

TEST_CASE("hash grid gets same impacts as brute force", "[hashgrid]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd displacements;
    vertices.resize(10, 2);
    vertices << -2.4375, 0.0, 0.0, 0.0, 2.4375, 0.0, -3.0, -3.0, 0.0, -3.0, 3.0,
        -3.0, -3.0, -6.0, -0.5625, -6.0, 0.5625, -6.0, 3.0, -6.0;

    edges.resize(9, 2);
    edges << 0, 1, 1, 2, 1, 4, 3, 4, 4, 5, 3, 6, 5, 9, 6, 7, 8, 9;

    displacements.resize(vertices.rows(), vertices.cols());
    for (int i = 0; i < 10; i++) {
        displacements.setRandom();
        displacements *= 3;

        EdgeVertexImpacts brute_force_ev_impacts;
        detect_edge_vertex_collisions(vertices, displacements, edges,
            brute_force_ev_impacts, DetectionMethod::BRUTE_FORCE);
        EdgeVertexImpacts hash_ev_impacts;
        detect_edge_vertex_collisions(vertices, displacements, edges,
            hash_ev_impacts, DetectionMethod::HASH_GRID);

        REQUIRE(brute_force_ev_impacts.size() == hash_ev_impacts.size());
        for (const auto& impact : brute_force_ev_impacts) {
            CHECK(std::find(
                      hash_ev_impacts.begin(), hash_ev_impacts.end(), impact)
                != hash_ev_impacts.end());
        }
    }
}
