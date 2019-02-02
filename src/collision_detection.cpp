// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <FixingCollisions/collision_detection.hpp>

namespace ccd {
ImpactPtr detect_edge_vertex_collision(
    const Eigen::VectorXd& /*vertex0_t0*/,
    const Eigen::VectorXd& /*vertex0_t1*/,
    const Eigen::VectorXd& /*edge_vertex1_t0*/,
    const Eigen::VectorXd& /*edge_vertex1_t1*/,
    const Eigen::VectorXd& /*edge_vertex2_t0*/,
    const Eigen::VectorXd& /*edge_vertex2_t1*/)
{
    throw "Not implemented yet";
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
    int dimensions = int(vertices_t0.cols());
    ImpactsPtr impacts(new Impacts());
    ImpactPtr potential_impact(nullptr);
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        auto edge = edges.block<1, 2>(edge_index, 0);
        for (int vertex_index = 0; vertex_index < vertices_t0.rows();
             vertex_index++) {
            if (vertex_index != edge(0) && vertex_index != edge(1)) {
                potential_impact = detect_edge_vertex_collision(
                    vertices_t0.block(vertex_index, 0, 1, dimensions),
                    vertices_t1.block(vertex_index, 0, 1, dimensions),
                    vertices_t0.block(edge(0), 0, 1, dimensions),
                    vertices_t1.block(edge(0), 0, 1, dimensions),
                    vertices_t0.block(edge(1), 0, 1, dimensions),
                    vertices_t1.block(edge(0), 0, 1, dimensions));
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
    throw "Hash Map collision detection is not implemented yet.";
}

}
