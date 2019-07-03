// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <ccd/collision_detection.hpp>

#include <autogen/time_of_impact_coeff.hpp>

#include <profiler.hpp>

namespace ccd {

// Determine if a single edge-vertext pair intersects.
void detect_edge_vertex_collisions_narrow_phase(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk, const EdgeVertexCandidate& ev_candidate,
    EdgeVertexImpacts& ev_impacts)
{
    double toi, alpha;
    PROFILE(
        bool are_colliding = compute_edge_vertex_time_of_impact(
            Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
        if (are_colliding) {
            ev_impacts.push_back(EdgeVertexImpact(toi, ev_candidate.edge_index,
                alpha, ev_candidate.vertex_index));
        },
        ProfiledPoint::DETECTING_COLLISIONS_NARROW_PHASE);
}

// Compute the time of impact of a point and edge moVing in 2D.
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk, double& toi, double& alpha,
    const double tolerance)
{
    // In 2D the intersecton of a point and edge moVing through time is the a
    // quadratic equation, at^2 + bt + c = 0. A sketch for the geometric proof
    // is: the point forms a line in space-time and the edge forms a bilinear
    // surface in space-time, so the impact is an intersection of the line and
    // bilinear surface. There can possibly 0, 1, 2, or infinite such
    // intersections. Therefore, the function for time of impact is a quadratic
    // equation. Importantly, this function only works for a 2D point and edge.
    // In order to extend this function to 3D one would have to solve a cubic
    // function. These coefficients were found using `python/intersections.py`
    // which uses sympy to compute a polynomial in terms of t.
    double a, b, c;
    autogen::time_of_impact_coeff<double>(Vi, Vj, Vk, Ui, Uj, Uk, a, b, c);

    auto check_solution = [&](double t, double& s) {
        return t >= 0 && t <= 1
            && temporal_parameterization_to_spatial(
                Vi, Vj, Vk, Ui, Uj, Uk, t, s, tolerance)
            && s >= 0 && s <= 1;
    };

    if (std::abs(a) > tolerance) { // Is the equation truly quadratic?
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
    } else if (std::abs(b) > tolerance) { // Is the equation truly linear?
        // Linear equation
        // bt + c = 0 => t = -c / b
        toi = -c / b;
        return check_solution(toi, alpha);
    } else if (std::abs(c) < tolerance) {
        // a = b = c = 0 => infinite solutions, but may not be on the edge.
        // Find the spatial locations along the line at t=0 and t=1
        double s0 = 0.0, s1 = 0.0;
        temporal_parameterization_to_spatial(
            Vi, Vj, Vk, Ui, Uj, Uk, 0, s0, tolerance);
        temporal_parameterization_to_spatial(
            Vi, Vj, Vk, Ui, Uj, Uk, 1, s1, tolerance);

        // Possible cases for trajectories:
        // - No impact to impact ():
        //      - s0 < 0 && s1 >= 0 -> sI = 0
        //      - s0 > 1 && s1 <= 1 -> sI = 1
        // - Impact to *: 0 <= s0 <= 1 -> sI = s0
        // - No impact to no impact: * -> NO_IMPACT
        if ((s0 < 0 && s1 >= 0) || (s0 > 1 && s1 <= 1)
            || (s0 >= 0 && s0 <= 1)) {
            alpha = s0 < 0 ? 0 : (s0 > 1 ? 1 : s0);
            toi = std::abs(s1 - s0) > tolerance ? ((alpha - s0) / (s1 - s0))
                                                : 0;
            return true;
        }
    }
    // else (c != 0) = 0 => no solution exists
    return false;
}

// Convert a temporal parameterization to a spatial parameterization.
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk, const double t, double& alpha,
    const double tolerance)
{
    Eigen::Vector2d numerator = Vi - Vk + t * (Ui - Uk);
    Eigen::Vector2d denominator = Vi - Vj + t * (Ui - Uj);
    assert(numerator.size() == denominator.size());

    if (std::abs(denominator.x()) > tolerance) {
        alpha = numerator.x() / denominator.x();
        return true;
    } else if (std::abs(denominator.y()) > tolerance) {
        alpha = numerator.y() / denominator.y();
        return true;
    } else if (numerator.isZero(tolerance)) {
        // The points are all equal at a time t.
        // I can prove n/d = 0/0 <=> p0(t) = p1(t) = p2(t).
        alpha = 0.5; // Any alpha will work, so I arbitrarily choose 0.5.
        return true;
    }
    return false;
}

} // namespace ccd
