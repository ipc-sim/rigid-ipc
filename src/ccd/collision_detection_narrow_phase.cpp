// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <cmath>

// Etienne Vouga's CCD using a root finder in floating points
#include <CTCD.h>

#include <ccd/time_of_impact.hpp>
#include <geometry/barycentric_coordinates.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Candidates& candidates,
    ConcurrentImpacts& impacts)
{
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    auto detect_ev_collision = [&](const EdgeVertexCandidate& ev_candidate) {
        long ei = ev_candidate.edge_index, vi = ev_candidate.vertex_index;
        double toi, alpha;
        bool is_impacting = detect_edge_vertex_collisions_narrow_phase(
            // Edge at t=0
            vertices.row(edges(ei, 0)), vertices.row(edges(ei, 1)),
            // Vertex at t=0
            vertices.row(vi),
            // Edge displacement (TODO: Convert this to edge at t=1)
            displacements.row(edges(ei, 0)), displacements.row(edges(ei, 1)),
            // Vertex displacement (TODO: Convert this to vertex at t=1)
            displacements.row(vi),
            // Output parameters
            toi, alpha);
        // Add the impact if the candidate is impacting
        if (is_impacting) {
            impacts.ev_impacts.emplace_back(toi, ei, alpha, vi);
        }
    };

    auto detect_ee_collision = [&](const EdgeEdgeCandidate& ee_candidate) {
        long e0i = ee_candidate.edge0_index, e1i = ee_candidate.edge1_index;
        double toi, alpha0, alpha1;
        bool is_impacting = detect_edge_edge_collisions_narrow_phase(
            // Edge 0 at t=0
            vertices.row(edges(e0i, 0)), vertices.row(edges(e0i, 1)),
            // Edge 1 at t=0
            vertices.row(edges(e1i, 0)), vertices.row(edges(e1i, 1)),
            // Edge 0 displacement (TODO: Convert this to edge at t=1)
            displacements.row(edges(e0i, 0)), displacements.row(edges(e0i, 1)),
            // Edge 1 displacement (TODO: Convert this to edge at t=1)
            displacements.row(edges(e1i, 0)), displacements.row(edges(e1i, 1)),
            // Output parameters
            toi, alpha0, alpha1);
        // Add the impact if the candidate is impacting
        if (is_impacting) {
            impacts.ee_impacts.emplace_back(toi, e0i, alpha0, e1i, alpha1);
        }
    };

    auto detect_fv_collision = [&](const FaceVertexCandidate& fv_candidate) {
        long fi = fv_candidate.face_index, vi = fv_candidate.vertex_index;
        double toi, u, v;
        bool is_impacting = detect_face_vertex_collisions_narrow_phase(
            // Face at t = 1
            vertices.row(faces(fi, 0)), vertices.row(faces(fi, 1)),
            vertices.row(faces(fi, 2)),
            // Vertex at t = 1
            vertices.row(fv_candidate.vertex_index),
            // Face displacement
            displacements.row(faces(fi, 0)), displacements.row(faces(fi, 1)),
            displacements.row(faces(fi, 2)),
            // Vertex displacement
            displacements.row(fv_candidate.vertex_index),
            // Output parameters
            toi, u, v);
        if (is_impacting) {
            impacts.fv_impacts.emplace_back(toi, fi, u, v, vi);
        }
    };

    tbb::parallel_invoke(
        [&] {
            impacts.ev_impacts.clear();
            tbb::parallel_for_each(
                candidates.ev_candidates, detect_ev_collision);
        },

        [&] {
            impacts.ee_impacts.clear();
            tbb::parallel_for_each(
                candidates.ee_candidates, detect_ee_collision);
        },

        [&] {
            impacts.fv_impacts.clear();
            tbb::parallel_for_each(
                candidates.fv_candidates, detect_fv_collision);
        });

    PROFILE_END(NARROW_PHASE);
    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    double& toi,
    double& alpha)
{
    return ccd::autodiff::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
}

bool detect_edge_edge_collisions_narrow_phase(
    const Eigen::Vector3d& Vi,
    const Eigen::Vector3d& Vj,
    const Eigen::Vector3d& Vk,
    const Eigen::Vector3d& Vl,
    const Eigen::Vector3d& Ui,
    const Eigen::Vector3d& Uj,
    const Eigen::Vector3d& Uk,
    const Eigen::Vector3d& Ul,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha)
{
    // TODO: Check if edges are parallel
    bool is_colliding = CTCD::edgeEdgeCTCD(
        Vi, Vj, Vk, Vl, Vi + Ui, Vj + Uj, Vk + Uk, Vl + Ul, /*eta=*/0, toi);
    if (is_colliding) {
        edge0_alpha = -1; // TODO: Compute this correctly
        edge1_alpha = -1; // TODO: Compute this correctly
    }
    return is_colliding;
}

bool detect_face_vertex_collisions_narrow_phase(
    const Eigen::Vector3d& Vi,
    const Eigen::Vector3d& Vj,
    const Eigen::Vector3d& Vk,
    const Eigen::Vector3d& Vl,
    const Eigen::Vector3d& Ui,
    const Eigen::Vector3d& Uj,
    const Eigen::Vector3d& Uk,
    const Eigen::Vector3d& Ul,
    double& toi,
    double& u,
    double& v)
{
    bool is_colliding = CTCD::vertexFaceCTCD(
        Vl, Vi, Vj, Vk, Vl + Ul, Vi + Ui, Vj + Uj, Vk + Uk, /*eta=*/0, toi);
    if (is_colliding) {
        // TODO: Consider moving this computation to an as needed basis
        double w;
        geometry::barycentric_coordinates(
            (Vl + Ul * toi).eval(), // vertex position at time of impact
            // first face vertex's position at time of impact
            (Vi + Ui * toi).eval(),
            // second face vertex's position at time of impact
            (Vj + Uj * toi).eval(),
            // third face vertex's position at time of impact
            (Vk + Uk * toi).eval(), u, v, w);
        assert(u + v + w == 1);
    }
    return is_colliding;
}

} // namespace ccd
