// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <FixingCollisions/collision_detection.hpp>
#include <FixingCollisions/degenerate_edge_error.hpp>
#include <FixingCollisions/not_implemented_error.hpp>

#include <iostream>

#define EPSILON (1e-8)

#define EDGE_VERTEX_PARAMS                                                     \
    vertex0, displacement0, edge_vertex1, edge_displacement1, edge_vertex2,    \
        edge_displacement2

namespace ccd {

// Convert a temporal parameterization to a spatial parameterization.
double temporal_parameterization_to_spatial(const Eigen::VectorXd& vertex0,
    const Eigen::VectorXd& displacement0, const Eigen::VectorXd& edge_vertex1,
    const Eigen::VectorXd& edge_displacement1,
    const Eigen::VectorXd& edge_vertex2,
    const Eigen::VectorXd& edge_displacement2, double t)
{
    Eigen::MatrixXd numerator
        = edge_vertex1 - vertex0 + t * (edge_displacement1 - displacement0);
    Eigen::MatrixXd denominator = edge_vertex1 - edge_vertex2
        + t * (edge_displacement1 - edge_displacement2);
    assert(numerator.size() == denominator.size());
    if (std::abs(denominator(0)) > EPSILON) {
        return numerator(0) / denominator(0);
    } else if (std::abs(denominator(1)) > EPSILON) {
        return numerator(1) / denominator(1);
    } else if (denominator.size() > 2 && std::abs(denominator(2)) > EPSILON) {
        return numerator(2) / denominator(2);
    }
    throw DegenerateEdgeError();
}

// Compute the time of impact of a point and edge moving in 2D.
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& vertex0,
    const Eigen::Vector2d& displacement0, const Eigen::Vector2d& edge_vertex1,
    const Eigen::Vector2d& edge_displacement1,
    const Eigen::Vector2d& edge_vertex2,
    const Eigen::Vector2d& edge_displacement2, double& toi, double& alpha)
{
    // In 2D the intersecton of a point and edge moving through time is the a
    // quadratic equation, at^2 + bt + c = 0. A sketch for the geometric proof
    // is: the point forms a line in space-time and the edge forms a
    // bilinear surface in space-time, so the impact is an intersection of the
    // line and bilinear surface. There can possibly 0, 1, 2, or infinite such
    // intersections. Therefore, the function for time of impact is a quadratic
    // equation.
    // Importantly, this function only works for a 2D point and edge. In order
    // to extend this function to 3D one would have to solve a cubic function.
    // These coefficients were found using `python/intersections.py` which uses
    // sympy to compute a polynomial in terms of t.
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

    if (std::abs(a) > EPSILON) { // Is the equation truly quadratic?
        // Quadratic equation
        // at^2 + bt + c = 0 => t = (-b ± sqrt(b^2 - 4ac)) / 2a
        double radicand = b * b - 4 * a * c;
        if (radicand >= 0) {
            double sqrt_rad = sqrt(radicand);
            // We know the time of impacts will be sorted earliest to latest.
            for (double sign : { -1, 1 }) {
                toi = (-b + sign * sqrt_rad) / (2 * a);
                try {
                    alpha = temporal_parameterization_to_spatial(
                        EDGE_VERTEX_PARAMS, toi);
                } catch (DegenerateEdgeError err) {
                    continue;
                }
                if (toi >= 0 && toi <= 1 && alpha >= 0 && alpha <= 1) {
                    return true;
                }
            }
        }
    } else if (std::abs(b) > EPSILON) { // Is the equation truly linear?
        // Linear equation
        // bt + c = 0 => t = -c / b
        toi = -c / b;
        alpha = temporal_parameterization_to_spatial(EDGE_VERTEX_PARAMS, toi);
        return toi >= 0 && toi <= 1 && alpha >= 0 && alpha <= 1;
    } else if (std::abs(c) < EPSILON) {
        // a = b = c = 0 => infinite solutions, but may not be on the edge.
        // Find the spatial locations along the line at t=0 and t=1
        double s0 = temporal_parameterization_to_spatial(EDGE_VERTEX_PARAMS, 0);
        double s1 = temporal_parameterization_to_spatial(EDGE_VERTEX_PARAMS, 1);

        // Possible cases for trajectories:
        // - No impact to impact ():
        //  - s0 < 0 && s1 >= 0 -> sI = 0
        //  - s0 > 1 && s1 <= 1 -> sI = 1
        // - impact to *: 0 <= s0 <= 1 -> sI = s0
        // - No impact to no impact: * -> NO_IMPACT
        if ((s0 < 0 && s1 >= 0) || (s0 > 1 && s1 <= 1)
            || (s0 >= 0 && s0 <= 1)) {
            alpha = s0 < 0 ? 0 : (s0 > 1 ? 1 : s0);
            toi = std::abs(s1 - s0) > EPSILON ? ((alpha - s0) / (s1 - s0)) : 0;
            return true;
        }
    }
    // else (c != 0) = 0 => no solution exists
    return false;
}

EdgeVertexImpactsPtr detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges, DetectionMethod method)
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

EdgeVertexImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges)
{
    EdgeVertexImpactsPtr impacts = std::make_shared<EdgeVertexImpacts>();
    double toi, alpha;
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        auto edge = edges.row(edge_index);
        for (int vertex_index = 0; vertex_index < vertices.rows();
             vertex_index++) {
            if (vertex_index != edge(0) && vertex_index != edge(1)) {
                if (compute_edge_vertex_time_of_impact(
                        vertices.row(vertex_index),
                        displacements.row(vertex_index), vertices.row(edge(0)),
                        displacements.row(edge(0)), vertices.row(edge(1)),
                        displacements.row(edge(1)), toi, alpha)) {
                    impacts->push_back(std::make_shared<EdgeVertexImpact>(
                        toi, edge_index, alpha, vertex_index));
                }
            }
        }
    }
    return impacts;
}

EdgeVertexImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& /* vertices */,
    const Eigen::MatrixXd& /* displacements */,
    const Eigen::MatrixX2i& /* edges */)
{
    throw NotImplementedError(
        "Hash Map collision detection is not implemented yet.");
}

}
