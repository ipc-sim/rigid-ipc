// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <FixingCollisions/collision_detection.hpp>

#include <iostream>

#include <FixingCollisions/not_implemented_error.hpp>

namespace ccd {
ImpactPtr detect_edge_vertex_collision(
    const Eigen::MatrixXd& /*vertex0_t0*/,
    const Eigen::MatrixXd& /*vertex0_t1*/,
    const Eigen::MatrixXd& /*edge_vertex1_t0*/,
    const Eigen::MatrixXd& /*edge_vertex1_t1*/,
    const Eigen::MatrixXd& /*edge_vertex2_t0*/,
    const Eigen::MatrixXd& /*edge_vertex2_t1*/)
{
    throw NotImplementedError("Detect edge-vertex collision not implemented");
}

ImpactsPtr detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::Matrix<int, Eigen::Dynamic, 2>& edges,
    DetectionMethod method)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
    switch (method) {
    case BRUTE_FORCE:
        return detect_edge_vertex_collisions_brute_force(
            vertices_t0, vertices_t1, edges);
    case HASH_MAP:
        return detect_edge_vertex_collisions_hash_map(
            vertices_t0, vertices_t1, edges);
    }
}

ImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::Matrix<int, Eigen::Dynamic, 2>& edges)
{
//    int dimensions = int(vertices_t0.cols());

    ImpactsPtr impacts(new Impacts());
    ImpactPtr potential_impact(nullptr);
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        auto edge = edges.middleRows<1>(edge_index);
        for (int vertex_index = 0; vertex_index < vertices_t0.rows();
             vertex_index++) {
            if (vertex_index != edge(0) && vertex_index != edge(1)) {
                potential_impact = detect_edge_vertex_collision(
                    vertices_t0.middleRows(vertex_index, 1),
                    vertices_t1.middleRows(vertex_index, 1),
                    vertices_t0.middleRows(edge(0), 1),
                    vertices_t1.middleRows(edge(0), 1),
                    vertices_t0.middleRows(edge(1), 1),
                    vertices_t1.middleRows(edge(1), 1));
                if (potential_impact) {
                    impacts->push_back(potential_impact);
                }
            }
        }
    }
    return impacts;
}

ImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& /*vertices_t0*/,
    const Eigen::MatrixXd& /*vertices_t1*/,
    const Eigen::Matrix<int, Eigen::Dynamic, 2>& /*edges*/)
{
    throw NotImplementedError("Hash Map collision detection is not implemented yet.");
}

void test()
{
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices_t0;
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> vertices_t1;
    Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
    vertices_t0 << 0, 2,
        1, 1,
        1, 3;
    vertices_t1 << 1, 2,
        0, 1,
        0, 3;
    edges << 1, 2;

    detect_edge_vertex_collisions(vertices_t0, vertices_t1, edges);
}
}
