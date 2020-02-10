#include <iostream>

#include <catch2/catch.hpp>

#include <ccd/collision_detection.hpp>
#include <ccd/time_of_impact.hpp>

using namespace ccd;

TEST_CASE("Test Continous Collision Detection", "[ccd][thisone]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    edges.row(0) << 1, 2;
    EdgeVertexImpacts ev_impacts;
    SECTION("Edge becomes degenerate")
    {
        vertices.row(0) << 0, 1;
        vertices.row(1) << -1, 0;
        vertices.row(2) << 1, 0;
        SECTION("Edge degenerates before impact")
        {
            displacements.row(0) << 0, -2;
            displacements.row(1) << 4, 0;
            displacements.row(2) << -4, 0;
            // The edge will become degenerate at t=0.25
        }
        // SECTION("Edge degenerates at impact")
        // {
        //     displacements.row(0) << 0, -2;
        //     displacements.row(1) << 2, 0;
        //     displacements.row(2) << -2, 0;
        // }
        SECTION("Edge degenerates after impact")
        {
            displacements.row(0) << 0, -2;
            displacements.row(1) << 1.5, 0;
            displacements.row(2) << -1.5, 0;
            // The edge will become degenerate at t=2/3
        }
        // The point will collide with the edge at t=0.5
        ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges, ev_impacts);
        CHECK(ev_impacts.size() == 1);
        if (ev_impacts.size() == 1) {
            CHECK(ev_impacts.front().time == Approx(0.5));
            CHECK(ev_impacts.front().alpha == Approx(0.5));
            CHECK(ev_impacts.front().vertex_index == 0);
            CHECK(ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Edge moving right; point moving left")
    {
        vertices.row(0) << -1, 0;
        vertices.row(1) << 1, -1;
        vertices.row(2) << 1, 1;
        displacements.row(0) << 2, 0;
        displacements.row(1) << -2, 0;
        displacements.row(2) << -2, 0;
        ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges, ev_impacts);
        CHECK(ev_impacts.size() == 1);
        if (ev_impacts.size() == 1) {
            CHECK(ev_impacts.front().time == Approx(0.5));
            CHECK(ev_impacts.front().alpha == Approx(0.5));
            CHECK(ev_impacts.front().vertex_index == 0);
            CHECK(ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Point on edge's line moving towards edge")
    {
        vertices.row(0) << 0, 0;
        vertices.row(1) << 0, 1;
        vertices.row(2) << 0, 2;
        displacements.row(0) << 0, 2;
        displacements.row(1) << 0, 0;
        displacements.row(2) << 0, 0;
        ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges, ev_impacts);
        CHECK(ev_impacts.size() == 1);
        if (ev_impacts.size() == 1) {
            CHECK(ev_impacts.front().time == Approx(0.5));
            CHECK(ev_impacts.front().alpha == Approx(0.0));
            CHECK(ev_impacts.front().vertex_index == 0);
            CHECK(ev_impacts.front().edge_index == 0);
        }
    }
    SECTION("Point and edge moving parallel")
    {
        vertices.row(0) << 0, 1;
        vertices.row(1) << 1, 0;
        vertices.row(2) << 1, 2;
        displacements.row(0) << 0, 1;
        displacements.row(1) << 0, 1;
        displacements.row(2) << 0, 1;
        ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges, ev_impacts);
        CHECK(ev_impacts.size() == 0);
    }
    SECTION("Point moving right; edge stretching vertically")
    {
        vertices.row(2) << 1, -1;
        displacements.row(2) << 0, -1;
        SECTION("Swap vertices order e_0 = [0, 2]")
        {
            vertices.row(0) << 1, 1;
            vertices.row(1) << 0, 0;
            displacements.row(0) << 0, 1;
            displacements.row(1) << 1, 0;
            edges(0) = 0;
            ccd::detect_edge_vertex_collisions(
                vertices, displacements, edges, ev_impacts);
            CHECK(ev_impacts.size() == 1);
            if (ev_impacts.size() == 1) {
                CHECK(ev_impacts.front().vertex_index == 1);
            }
        }
        SECTION("Swap vertices order e_0 = [1, 2]")
        {
            vertices.row(0) << 0, 0;
            vertices.row(1) << 1, 1;
            displacements.row(0) << 1, 0;
            displacements.row(1) << 0, 1;
            edges(0) = 1;
            ccd::detect_edge_vertex_collisions(
                vertices, displacements, edges, ev_impacts);
            CHECK(ev_impacts.size() == 1);
            if (ev_impacts.size() == 1) {
                CHECK(ev_impacts.front().vertex_index == 0);
            }
        }
        CHECK(ev_impacts.front().edge_index == 0);
        CHECK(ev_impacts.front().time == Approx(1.0));
        CHECK(ev_impacts.front().alpha == Approx(0.5));
    }
}
