#include <iostream>

#include <catch2/catch.hpp>

#include <FixingCollisions/collision_detection.hpp>

using namespace ccd;

TEST_CASE("ToI test 1", "[ccd]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    vertices.row(0) << -1, 0;
    vertices.row(1) << 1, -1;
    vertices.row(2) << 1, 1;
    displacements.row(0) << 2, 0;
    displacements.row(1) << -2, 0;
    displacements.row(2) << -2, 0;
    edges.row(0) << 1, 2;
    ccd::ImpactsPtr impacts
        = ccd::detect_edge_vertex_collisions(vertices, displacements, edges);
    REQUIRE(impacts->size() == 1);
    REQUIRE(impacts->front()->time == Approx(0.5));
    REQUIRE(impacts->front()->alpha == Approx(0.5));
}

TEST_CASE("ToI test 2", "[ccd]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    vertices.row(0) << 0, 0;
    vertices.row(1) << 0, 1;
    vertices.row(2) << 0, 2;
    displacements.row(0) << 0, 2;
    displacements.row(1) << 0, 0;
    displacements.row(2) << 0, 0;
    edges.row(0) << 1, 2;
    ccd::ImpactsPtr impacts
        = ccd::detect_edge_vertex_collisions(vertices, displacements, edges);
    REQUIRE(impacts->size() == 1);
    REQUIRE(impacts->front()->time == Approx(0.5));
    REQUIRE(impacts->front()->alpha == Approx(0.0));
}

TEST_CASE("ToI test 3", "[ccd]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    vertices.row(0) << 0, 1;
    vertices.row(1) << 1, 0;
    vertices.row(2) << 1, 2;
    displacements.row(0) << 0, 1;
    displacements.row(1) << 0, 1;
    displacements.row(2) << 0, 1;
    edges.row(0) << 1, 2;
    ccd::ImpactsPtr impacts
        = ccd::detect_edge_vertex_collisions(vertices, displacements, edges);
    REQUIRE(impacts->size() == 0);
}

TEST_CASE("ToI test 4", "[ccd]")
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> displacements;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    vertices.row(0) << 0, 0;
    vertices.row(1) << 1, 1;
    vertices.row(2) << 1, -1;
    displacements.row(0) << 1, 0;
    displacements.row(1) << 0, 1;
    displacements.row(2) << 0, -1;
    edges.row(0) << 1, 2;
    ccd::ImpactsPtr impacts
        = ccd::detect_edge_vertex_collisions(vertices, displacements, edges);
    REQUIRE(impacts->size() == 1);
    REQUIRE(impacts->front()->time == Approx(1.0));
    REQUIRE(impacts->front()->alpha == Approx(0.5));
}
