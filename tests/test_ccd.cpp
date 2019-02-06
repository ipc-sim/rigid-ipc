#include <iostream>

#include <FixingCollisions/collision_detection.hpp>

using namespace ccd;

void test1()
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
    for (auto const& impact : *impacts) {
        std::cout << "Impact between point " << impact->vertex_index
                  << " and edge " << impact->edge_index << " at time "
                  << impact->time << std::endl;
    }
}

void test2()
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
    for (auto const& impact : *impacts) {
        std::cout << "Impact between point " << impact->vertex_index
                  << " and edge " << impact->edge_index << " at time "
                  << impact->time << std::endl;
    }
}

void test3()
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
    for (auto const& impact : *impacts) {
        std::cout << "Impact between point " << impact->vertex_index
                  << " and edge " << impact->edge_index << " at time "
                  << impact->time << std::endl;
    }
}

void test4()
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
    for (auto const& impact : *impacts) {
        std::cout << "Impact between point " << impact->vertex_index
                  << " and edge " << impact->edge_index << " at time "
                  << impact->time << std::endl;
    }
}
