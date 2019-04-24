// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <ccd/collision_detection.hpp>
#include <ccd/degenerate_edge_error.hpp>
#include <ccd/not_implemented_error.hpp>

#include <iostream>

#define EPSILON (1e-8)

namespace ccd {

// Convert a temporal parameterization to a spatial parameterization.
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& vertex0,
    const Eigen::Vector2d& displacement0, const Eigen::Vector2d& edge_vertex1,
    const Eigen::Vector2d& edge_displacement1,
    const Eigen::Vector2d& edge_vertex2,
    const Eigen::Vector2d& edge_displacement2, const double t, double& alpha)
{
    Eigen::Vector2d numerator
        = edge_vertex1 - vertex0 + t * (edge_displacement1 - displacement0);
    Eigen::Vector2d denominator = edge_vertex1 - edge_vertex2
        + t * (edge_displacement1 - edge_displacement2);
    assert(numerator.size() == denominator.size());

    if (std::abs(denominator(0)) > EPSILON) {
        alpha = numerator(0) / denominator(0);
        return true;
    } else if (std::abs(denominator(1)) > EPSILON) {
        alpha = numerator(1) / denominator(1);
        return true;
    } else if (numerator.isZero(EPSILON)) {
        // The points are all equal at a time t.
        // I can prove n/d = 0/0 <=> p0(t) = p1(t) = p2(t).
        alpha = 0.5; // Any alpha will work, so I arbitrarily choose 0.5.
        return true;
    }
    return false;
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
    // is: the point forms a line in space-time and the edge forms a bilinear
    // surface in space-time, so the impact is an intersection of the line and
    // bilinear surface. There can possibly 0, 1, 2, or infinite such
    // intersections. Therefore, the function for time of impact is a quadratic
    // equation. Importantly, this function only works for a 2D point and edge.
    // In order to extend this function to 3D one would have to solve a cubic
    // function. These coefficients were found using `python/intersections.py`
    // which uses sympy to compute a polynomial in terms of t.
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

    auto check_solution = [&](double t, double& s) {
        return t >= 0 && t <= 1
            && temporal_parameterization_to_spatial(vertex0, displacement0,
                edge_vertex1, edge_displacement1, edge_vertex2,
                edge_displacement2, t, s)
            && s >= 0 && s <= 1;
    };

    if (std::abs(a) > EPSILON) { // Is the equation truly quadratic?
        // Quadratic equation
        // at^2 + bt + c = 0 => t = (-b ± sqrt(b^2 - 4ac)) / 2a
        double radicand = b * b - 4 * a * c;
        if (radicand >= 0) {
            double sqrt_rad = sqrt(radicand);
            double toi0 = (-b + sqrt_rad) / (2 * a);
            double toi1 = (-b - sqrt_rad) / (2 * a);
            double alpha0, alpha1;
            bool is_toi0_valid = check_solution(toi0, alpha0);
            bool is_toi1_valid = check_solution(toi1, alpha1);

            if (is_toi0_valid) {
                toi = (!is_toi1_valid || toi0 < toi1) ? toi0 : toi1;
                alpha = (!is_toi1_valid || toi0 < toi1) ? alpha0 : alpha1;
                return true;
            } else if (is_toi1_valid) {
                toi = toi1;
                alpha = alpha1;
                return true;
            }
        }
    } else if (std::abs(b) > EPSILON) { // Is the equation truly linear?
        // Linear equation
        // bt + c = 0 => t = -c / b
        toi = -c / b;
        return check_solution(toi, alpha);
    } else if (std::abs(c) < EPSILON) {
        // a = b = c = 0 => infinite solutions, but may not be on the edge.
        // Find the spatial locations along the line at t=0 and t=1
        double s0(0.0), s1(0.0);
        temporal_parameterization_to_spatial(vertex0, displacement0,
            edge_vertex1, edge_displacement1, edge_vertex2, edge_displacement2,
            0, s0);
        temporal_parameterization_to_spatial(vertex0, displacement0,
            edge_vertex1, edge_displacement1, edge_vertex2, edge_displacement2,
            1, s1);

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

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts, DetectionMethod method, bool reset_impacts)
{
    assert(vertices.size() == displacements.size());
    assert(method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_MAP);
    switch (method) {
    case BRUTE_FORCE:
        detect_edge_vertex_collisions_brute_force(
            vertices, displacements, edges, ev_impacts, reset_impacts);
        break;
    case HASH_MAP:
        detect_edge_vertex_collisions_hash_map(
            vertices, displacements, edges, ev_impacts, reset_impacts);
        break;
    }
}

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_edge_vertex_collisions_brute_force(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts, bool reset_impacts)
{
    double toi = -1, alpha = -1;

    if (reset_impacts) {
        ev_impacts.clear();
    }

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        was_impact_found = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::RowMajor>::Constant(edges.rows(), vertices.rows(), false);

    // If we do not recompute vertex impacts then we need to prevent duplicates
    for (EdgeVertexImpact ev_impact : ev_impacts) {
        was_impact_found(ev_impact.edge_index, ev_impact.vertex_index) = true;
    }

    EdgeVertexImpact ev_impact;
    // Loop over all edges
    for (int edge_idx = 0; edge_idx < edges.rows(); edge_idx++) {
        // Get the edge endpoint indices
        Eigen::Vector2i edge = edges.row(edge_idx);

        // Loop over all vertices
        for (int vertex_idx = 0; vertex_idx < vertices.rows(); vertex_idx++) {
            // Check that the vertex is not an endpoint of the edge
            if (vertex_idx != edge(0) && vertex_idx != edge(1)) {
                // Check if there is a collision between the vertex and edge
                if (!was_impact_found(edge_idx, vertex_idx)
                    && compute_edge_vertex_time_of_impact(
                        vertices.row(vertex_idx), displacements.row(vertex_idx),
                        vertices.row(edge(0)), displacements.row(edge(0)),
                        vertices.row(edge(1)), displacements.row(edge(1)), toi,
                        alpha)) {

                    ev_impact.time = toi;
                    ev_impact.edge_index = edge_idx;
                    ev_impact.alpha = alpha;
                    ev_impact.vertex_index = vertex_idx;

                    ev_impacts.push_back(ev_impact);
                }
            }
        }
    }
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& /* vertices */,
    const Eigen::MatrixXd& /* displacements */,
    const Eigen::MatrixX2i& /* edges */, EdgeVertexImpacts& /* ev_impacts */,
    bool /* reset_impacts */)
{
    throw NotImplementedError(
        "Hash Map collision detection is not implemented yet.");
}

} // namespace ccd
