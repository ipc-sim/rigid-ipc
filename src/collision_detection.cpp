// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <FixingCollisions/collision_detection.hpp>

#include <iostream>

#include <FixingCollisions/not_implemented_error.hpp>

#define NO_IMPACT (-1)
#define EPSILON (1e-8)

namespace ccd {

double temporal_parameterization_to_spactial(const Eigen::MatrixXd& vertex0,
    const Eigen::MatrixXd& displacement0, const Eigen::MatrixXd& edge_vertex1,
    const Eigen::MatrixXd& edge_displacement1,
    const Eigen::MatrixXd& edge_vertex2,
    const Eigen::MatrixXd& edge_displacement2, double t)
{
    Eigen::MatrixXd numerator
        = edge_vertex1 - vertex0 + t * (edge_displacement1 - displacement0);
    Eigen::MatrixXd denominator = edge_vertex1 - edge_vertex2
        + t * (edge_displacement1 - edge_displacement2);
    if (std::abs(denominator(0)) > EPSILON) {
        return numerator(0) / denominator(0);
    } else if (std::abs(denominator(1)) > EPSILON) {
        return numerator(1) / denominator(1);
    }
    throw "Edge is a singular non-moving point.";
}

double compute_edge_vertex_time_of_impact(const Eigen::MatrixX2d& vertex0,
    const Eigen::MatrixX2d& displacement0, const Eigen::MatrixX2d& edge_vertex1,
    const Eigen::MatrixX2d& edge_displacement1,
    const Eigen::MatrixX2d& edge_vertex2,
    const Eigen::MatrixX2d& edge_displacement2)
{
    // a*t^2 + b*t + c = 0
    double a
        = displacement0(0) * (edge_displacement2(1) - edge_displacement1(1))
        + displacement0(1) * (edge_displacement1(0) - edge_displacement2(0))
        - edge_displacement1(0) * edge_displacement2(1)
        + edge_displacement1(1) * edge_displacement2(0);
    double b = vertex0(0) * (edge_displacement2(1) - edge_displacement1(1))
        + vertex0(1) * (edge_displacement1(0) - edge_displacement2(0))
        + edge_vertex1(0) * (displacement0(1) - edge_displacement2(1))
        + edge_vertex1(1) * (edge_displacement2(0) - displacement0(0))
        + edge_vertex2(0) * (edge_displacement1(1) - displacement0(1))
        + edge_vertex2(1) * (displacement0(0) - edge_displacement1(0));
    double c = vertex0(0) * (edge_vertex2(1) - edge_vertex1(1))
        + vertex0(1) * (edge_vertex1(0) - edge_vertex2(0))
        + edge_vertex1(1) * edge_vertex2(0) - edge_vertex1(0) * edge_vertex2(1);

    if (std::abs(a) > EPSILON) {
        // Quadratic equation
        double radicand = b * b - 4 * a * c;
        if (radicand >= 0) {
            // TODO: Check that a is not zero
            double tmp0 = -b / (2 * a);
            double tmp1 = sqrt(radicand) / (2 * a);
            double t0 = tmp0 + tmp1;
            double t1 = tmp0 - tmp1;
            double s0 = temporal_parameterization_to_spactial(vertex0,
                displacement0, edge_vertex1, edge_displacement1, edge_vertex2,
                edge_displacement2, t0);
            double s1 = temporal_parameterization_to_spactial(vertex0,
                displacement0, edge_vertex1, edge_displacement1, edge_vertex2,
                edge_displacement2, t1);
            bool is_t0_valid = t0 >= 0 && t0 <= 1 && s0 >= 0 && s0 <= 1;
            bool is_t1_valid = t1 >= 0 && t1 <= 1 && s1 >= 0 && s1 <= 1;
            if (is_t0_valid) {
                return is_t1_valid ? std::min(t0, t1) : t0;
            } else {
                return is_t1_valid ? t1 : NO_IMPACT;
            }
        }
        return NO_IMPACT;
    } else if (std::abs(b) > EPSILON) {
        // Linear equation
        // b * t + c = 0 => t = -c / b
        double t = -c / b;
        double s = temporal_parameterization_to_spactial(vertex0, displacement0,
            edge_vertex1, edge_displacement1, edge_vertex2, edge_displacement2,
            t);
        return (t >= 0 && t <= 1 && s >= 0 && s <= 1) ? t : NO_IMPACT;
    } else if (std::abs(c) > EPSILON) {
        // (c != 0) = 0 => no solution exists
        return NO_IMPACT;
    } else {
        // a = b = c = 0 => infinite solutions, but may not be on the edge.
        throw NotImplementedError("Case a = b = c = 0 is not implemented yet.");
    }
}

ImpactsPtr detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    DetectionMethod method)
{
    assert(vertices.size() == displacements.size());
    switch (method) {
    case BRUTE_FORCE:
        return detect_edge_vertex_collisions_brute_force(
            vertices, displacements, edges);
    case HASH_MAP:
        return detect_edge_vertex_collisions_hash_map(
            vertices, displacements, edges);
    }
}

ImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges)
{
    ImpactsPtr impacts(new Impacts());
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        auto edge = edges.row(edge_index);
        for (int vertex_index = 0; vertex_index < vertices.rows();
             vertex_index++) {
            if (vertex_index != edge(0) && vertex_index != edge(1)) {
                double potential_time_of_impact
                    = compute_edge_vertex_time_of_impact(
                        vertices.row(vertex_index),
                        displacements.row(vertex_index), vertices.row(edge(0)),
                        displacements.row(edge(0)), vertices.row(edge(1)),
                        displacements.row(edge(1)));
                if (potential_time_of_impact != NO_IMPACT) {
                    impacts->push_back(ImpactPtr(new Impact()));
                    impacts->back()->vertex_index = vertex_index;
                    impacts->back()->edge_index = edge_index;
                    impacts->back()->time = potential_time_of_impact;
                }
            }
        }
    }
    return impacts;
}

ImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& /* vertices */,
    const Eigen::MatrixXd& /* displacements */,
    const Eigen::MatrixX2i& /* edges */)
{
    throw NotImplementedError(
        "Hash Map collision detection is not implemented yet.");
}

void test()
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
}
