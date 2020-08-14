// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <cmath>

// Etienne Vouga's CCD using a root finder in floating points
#include <CTCD.h>

#include <ccd/time_of_impact.hpp>
#include <constants.hpp>
#include <geometry/barycentric_coordinates.hpp>
#include <geometry/projection.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Candidates& candidates,
    ConcurrentImpacts& impacts)
{
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);
    PROFILE_START(NARROW_PHASE);

    auto detect_ev_collision = [&](const EdgeVertexCandidate& ev_candidate) {
        long ei = ev_candidate.edge_index, vi = ev_candidate.vertex_index;
        double toi, alpha;
        bool is_impacting = detect_edge_vertex_collisions_narrow_phase(
            // Edge at t=0
            vertices_t0.row(edges(ei, 0)), vertices_t0.row(edges(ei, 1)),
            // Vertex at t=0
            vertices_t0.row(vi),
            // Edge at t=1
            vertices_t1.row(edges(ei, 0)), vertices_t1.row(edges(ei, 1)),
            // Vertex at t=1
            vertices_t1.row(vi),
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
            vertices_t0.row(edges(e0i, 0)), vertices_t0.row(edges(e0i, 1)),
            // Edge 1 at t=0
            vertices_t0.row(edges(e1i, 0)), vertices_t0.row(edges(e1i, 1)),
            // Edge 0 at t=1
            vertices_t1.row(edges(e0i, 0)), vertices_t1.row(edges(e0i, 1)),
            // Edge 1 at t=1
            vertices_t1.row(edges(e1i, 0)), vertices_t1.row(edges(e1i, 1)),
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
            // Face at t = 0
            vertices_t0.row(faces(fi, 0)), vertices_t0.row(faces(fi, 1)),
            vertices_t0.row(faces(fi, 2)),
            // Vertex at t = 0
            vertices_t0.row(fv_candidate.vertex_index),
            // Face at t = 1
            vertices_t1.row(faces(fi, 0)), vertices_t1.row(faces(fi, 1)),
            vertices_t1.row(faces(fi, 2)),
            // Vertex at t = 1
            vertices_t1.row(fv_candidate.vertex_index),
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
}

// Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const Eigen::Vector2d& edge_vertex0_t0,
    const Eigen::Vector2d& edge_vertex1_t0,
    const Eigen::Vector2d& vertex_t0,
    const Eigen::Vector2d& edge_vertex0_t1,
    const Eigen::Vector2d& edge_vertex1_t1,
    const Eigen::Vector2d& vertex_t1,
    double& toi,
    double& alpha)
{
    // Expects the arguments as position and displacements
    return ccd::autodiff::compute_edge_vertex_time_of_impact(
        edge_vertex0_t0, edge_vertex1_t0, vertex_t0,
        (edge_vertex0_t1 - edge_vertex0_t0).eval(),
        (edge_vertex1_t1 - edge_vertex1_t0).eval(),
        (vertex_t1 - vertex_t0).eval(), toi, alpha);
}

bool detect_edge_edge_collisions_narrow_phase(
    const Eigen::Vector3d& edge0_vertex0_t0,
    const Eigen::Vector3d& edge0_vertex1_t0,
    const Eigen::Vector3d& edge1_vertex0_t0,
    const Eigen::Vector3d& edge1_vertex1_t0,
    const Eigen::Vector3d& edge0_vertex0_t1,
    const Eigen::Vector3d& edge0_vertex1_t1,
    const Eigen::Vector3d& edge1_vertex0_t1,
    const Eigen::Vector3d& edge1_vertex1_t1,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha)
{
    // TODO: Check if edges are parallel
    bool is_colliding = CTCD::edgeEdgeCTCD(
        edge0_vertex0_t0, edge0_vertex1_t0, edge1_vertex0_t0, edge1_vertex1_t0,
        edge0_vertex0_t1, edge0_vertex1_t1, edge1_vertex0_t1, edge1_vertex1_t1,
        Constants::LINEARIZED_CCD_ETA, toi);
    if (is_colliding) {
        Eigen::Vector3d edge0_vertex0_toi =
            (edge0_vertex0_t1 - edge0_vertex0_t0) * toi + edge0_vertex0_t0;
        Eigen::Vector3d edge0_vertex1_toi =
            (edge0_vertex1_t1 - edge0_vertex1_t0) * toi + edge0_vertex1_t0;
        Eigen::Vector3d edge1_vertex0_toi =
            (edge1_vertex0_t1 - edge1_vertex0_t0) * toi + edge1_vertex0_t0;
        Eigen::Vector3d edge1_vertex1_toi =
            (edge1_vertex1_t1 - edge1_vertex1_t0) * toi + edge1_vertex1_t0;
        geometry::project_segment_to_segment(
            edge0_vertex0_toi, edge0_vertex1_toi, edge1_vertex0_toi,
            edge1_vertex1_toi, edge0_alpha, edge1_alpha);
        assert(
            edge0_alpha > -Constants::PARAMETER_ASSERTION_TOL
            && edge0_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
        assert(
            edge1_alpha > -Constants::PARAMETER_ASSERTION_TOL
            && edge1_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
    }
    return is_colliding;
}

bool detect_face_vertex_collisions_narrow_phase(
    const Eigen::Vector3d& face_vertex0_t0,
    const Eigen::Vector3d& face_vertex1_t0,
    const Eigen::Vector3d& face_vertex2_t0,
    const Eigen::Vector3d& vertex_t0,
    const Eigen::Vector3d& face_vertex0_t1,
    const Eigen::Vector3d& face_vertex1_t1,
    const Eigen::Vector3d& face_vertex2_t1,
    const Eigen::Vector3d& vertex_t1,
    double& toi,
    double& u,
    double& v)
{
    bool is_colliding = CTCD::vertexFaceCTCD(
        vertex_t0, face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, //
        vertex_t1, face_vertex0_t1, face_vertex1_t1, face_vertex2_t1,
        Constants::LINEARIZED_CCD_ETA, toi);
    if (is_colliding) {
        // TODO: Consider moving this computation to an as needed basis
        Eigen::Vector3d face_vertex0_toi =
            (face_vertex0_t1 - face_vertex0_t0) * toi + face_vertex0_t0;
        Eigen::Vector3d face_vertex1_toi =
            (face_vertex1_t1 - face_vertex1_t0) * toi + face_vertex1_t0;
        Eigen::Vector3d face_vertex2_toi =
            (face_vertex2_t1 - face_vertex2_t0) * toi + face_vertex2_t0;
        Eigen::Vector3d vertex1_toi = (vertex_t1 - vertex_t0) * toi + vertex_t0;
        double w;
        geometry::barycentric_coordinates(
            vertex1_toi, face_vertex0_toi, face_vertex1_toi, face_vertex2_toi,
            u, v, w);
        assert(u + v + w < 1 + 1e-12 && u + v + w > 1 - 1e-12);
    }
    return is_colliding;
}

} // namespace ccd
