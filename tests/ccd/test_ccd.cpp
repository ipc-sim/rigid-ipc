#include <catch2/catch.hpp>

#include <ccd/linear/ccd.hpp>
#include <ccd/linear/edge_vertex_ccd.hpp>

using namespace ccd;

TEST_CASE("Test Continous Collision Detection", "[ccd][thisone]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices_t0, vertices_t1;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    edges.row(0) << 1, 2;
    ConcurrentImpacts impacts;
    SECTION("Edge becomes degenerate")
    {
        vertices_t0.row(0) << 0, 1;
        vertices_t0.row(1) << -1, 0;
        vertices_t0.row(2) << 1, 0;
        SECTION("Edge degenerates before impact")
        {
            vertices_t1.row(0) << 0, -1;
            vertices_t1.row(1) << 3, 0;
            vertices_t1.row(2) << -3, 0;
            // The edge will become degenerate at t=0.25
        }
        SECTION("Edge degenerates after impact")
        {
            vertices_t1.row(0) << 0, -1;
            vertices_t1.row(1) << 0.5, 0;
            vertices_t1.row(2) << -0.5, 0;
            // The edge will become degenerate at t=2/3
        }
        // The point will collide with the edge at t=0.5
        ccd::detect_collisions(
            vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
        CHECK(impacts.ev_impacts.size() == 1);
        if (impacts.ev_impacts.size() == 1) {
            CHECK(impacts.ev_impacts.front().time == Approx(0.5));
            CHECK(impacts.ev_impacts.front().alpha == Approx(0.5));
            CHECK(impacts.ev_impacts.front().vertex_index == 0);
            CHECK(impacts.ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Edge moving right; point moving left")
    {
        vertices_t0.row(0) << -1, 0;
        vertices_t0.row(1) << 1, -1;
        vertices_t0.row(2) << 1, 1;

        vertices_t1.row(0) << 1, 0;
        vertices_t1.row(1) << -1, -1;
        vertices_t1.row(2) << -1, 1;

        ccd::detect_collisions(
            vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
        CHECK(impacts.ev_impacts.size() == 1);
        if (impacts.ev_impacts.size() == 1) {
            CHECK(impacts.ev_impacts.front().time == Approx(0.5));
            CHECK(impacts.ev_impacts.front().alpha == Approx(0.5));
            CHECK(impacts.ev_impacts.front().vertex_index == 0);
            CHECK(impacts.ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Point on edge's line moving towards edge")
    {
        vertices_t0.row(0) << 0, 0;
        vertices_t0.row(1) << 0, 1;
        vertices_t0.row(2) << 0, 2;

        vertices_t1.row(0) << 0, 2;
        vertices_t1.row(1) << 0, 1;
        vertices_t1.row(2) << 0, 2;

        ccd::detect_collisions(
            vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
        CHECK(impacts.ev_impacts.size() == 1);
        if (impacts.ev_impacts.size() == 1) {
            CHECK(impacts.ev_impacts.front().time == Approx(0.5));
            CHECK(impacts.ev_impacts.front().alpha == Approx(0.0));
            CHECK(impacts.ev_impacts.front().vertex_index == 0);
            CHECK(impacts.ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Point and edge moving parallel")
    {
        vertices_t0.row(0) << 0, 1;
        vertices_t0.row(1) << 1, 0;
        vertices_t0.row(2) << 1, 2;

        vertices_t1.row(0) << 0, 2;
        vertices_t1.row(1) << 1, 1;
        vertices_t1.row(2) << 1, 3;

        ccd::detect_collisions(
            vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
            Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
        CHECK(impacts.ev_impacts.size() == 0);
    }
    SECTION("Point moving right; edge stretching vertically")
    {
        vertices_t0.row(2) << 1, -1;
        vertices_t1.row(2) << 1, -2;
        SECTION("Swap vertices order e_0 = [0, 2]")
        {
            vertices_t0.row(0) << 1, 1;
            vertices_t0.row(1) << 0, 0;

            vertices_t1.row(0) << 1, 2;
            vertices_t1.row(1) << 1, 0;

            edges(0) = 0;
            ccd::detect_collisions(
                vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
                Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
            CHECK(impacts.ev_impacts.size() == 1);
            if (impacts.ev_impacts.size() == 1) {
                CHECK(impacts.ev_impacts.front().vertex_index == 1);
            }
        }
        SECTION("Swap vertices order e_0 = [1, 2]")
        {
            vertices_t0.row(0) << 0, 0;
            vertices_t0.row(1) << 1, 1;

            vertices_t1.row(0) << 1, 0;
            vertices_t1.row(1) << 1, 2;

            edges(0) = 1;
            ccd::detect_collisions(
                vertices_t0, vertices_t1, edges, Eigen::MatrixXi(0, 3),
                Eigen::VectorXi(), ccd::CollisionType::EDGE_VERTEX, impacts);
            CHECK(impacts.ev_impacts.size() == 1);
            if (impacts.ev_impacts.size() == 1) {
                CHECK(impacts.ev_impacts.front().vertex_index == 0);
            }
        }
        CHECK(impacts.ev_impacts.front().edge_index == 0);
        CHECK(impacts.ev_impacts.front().time == Approx(1.0));
        CHECK(impacts.ev_impacts.front().alpha == Approx(0.5));
    }
}
