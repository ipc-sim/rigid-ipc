#include <iostream>

#include <catch2/catch.hpp>

#include <FixingCollisions/collision_detection.hpp>

using namespace ccd;

TEST_CASE("Test Continous Collision Detection", "[ccd]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    edges.row(0) << 1, 2;
    SECTION("Edge moving right; point moving left")
    {
        vertices.row(0) << -1, 0;
        vertices.row(1) << 1, -1;
        vertices.row(2) << 1, 1;
        displacements.row(0) << 2, 0;
        displacements.row(1) << -2, 0;
        displacements.row(2) << -2, 0;
        auto impacts = ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges);
        REQUIRE(impacts->size() == 1);
        REQUIRE(impacts->front()->time == Approx(0.5));
        REQUIRE(impacts->front()->alpha == Approx(0.5));
        REQUIRE(impacts->front()->vertex_index == 0);
        REQUIRE(impacts->front()->edge_index == 0);
    }
    SECTION("Point on edge's line moving towards edge")
    {
        vertices.row(0) << 0, 0;
        vertices.row(1) << 0, 1;
        vertices.row(2) << 0, 2;
        displacements.row(0) << 0, 2;
        displacements.row(1) << 0, 0;
        displacements.row(2) << 0, 0;
        auto impacts = ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges);
        REQUIRE(impacts->size() == 1);
        REQUIRE(impacts->front()->time == Approx(0.5));
        REQUIRE(impacts->front()->alpha == Approx(0.0));
        REQUIRE(impacts->front()->vertex_index == 0);
        REQUIRE(impacts->front()->edge_index == 0);
    }
    SECTION("Point and edge moving parallel")
    {
        vertices.row(0) << 0, 1;
        vertices.row(1) << 1, 0;
        vertices.row(2) << 1, 2;
        displacements.row(0) << 0, 1;
        displacements.row(1) << 0, 1;
        displacements.row(2) << 0, 1;
        auto impacts = ccd::detect_edge_vertex_collisions(
            vertices, displacements, edges);
        REQUIRE(impacts->size() == 0);
    }
    SECTION("Point moving right; edge stretching vertically")
    {
        EdgeVertexImpactsPtr impacts;
        vertices.row(2) << 1, -1;
        displacements.row(2) << 0, -1;
        SECTION("Swap vertices order e_0 = [0, 2]")
        {
            vertices.row(0) << 1, 1;
            vertices.row(1) << 0, 0;
            displacements.row(0) << 0, 1;
            displacements.row(1) << 1, 0;
            edges(0) = 0;
            impacts = ccd::detect_edge_vertex_collisions(
                vertices, displacements, edges);
            REQUIRE(impacts->front()->vertex_index == 1);
        }
        SECTION("Swap vertices order e_0 = [1, 2]")
        {
            vertices.row(0) << 0, 0;
            vertices.row(1) << 1, 1;
            displacements.row(0) << 1, 0;
            displacements.row(1) << 0, 1;
            edges(0) = 1;
            impacts = ccd::detect_edge_vertex_collisions(
                vertices, displacements, edges);
            REQUIRE(impacts->front()->vertex_index == 0);
        }
        REQUIRE(impacts->front()->edge_index == 0);
        REQUIRE(impacts->size() == 1);
        REQUIRE(impacts->front()->time == Approx(1.0));
        REQUIRE(impacts->front()->alpha == Approx(0.5));
    }
}
