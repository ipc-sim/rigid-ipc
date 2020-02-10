// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

// Ettien Vouga's CCD using a root finder in floating points
#include <CTCD.h>
#include <ccd/time_of_impact.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const EdgeVertexCandidates& ev_candidates,
    const EdgeEdgeCandidates& ee_candidates,
    const FaceVertexCandidates& fv_candidates,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts)
{
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    tbb::parallel_invoke(
        [&] {
            ev_impacts.clear();
            tbb::parallel_for_each(
                ev_candidates, [&](const EdgeVertexCandidate& ev_candidate) {
                    detect_edge_vertex_collisions_narrow_phase(
                        vertices.row(edges(ev_candidate.edge_index, 0)),
                        vertices.row(edges(ev_candidate.edge_index, 1)),
                        vertices.row(ev_candidate.vertex_index),
                        displacements.row(edges(ev_candidate.edge_index, 0)),
                        displacements.row(edges(ev_candidate.edge_index, 1)),
                        displacements.row(ev_candidate.vertex_index),
                        ev_candidate, ev_impacts);
                });
        },

        [&] {
            ee_impacts.clear();
            tbb::parallel_for_each(
                ee_candidates, [&](const EdgeEdgeCandidate& ee_candidate) {
                    detect_edge_edge_collisions_narrow_phase(
                        vertices.row(edges(ee_candidate.edge0_index, 0)),
                        vertices.row(edges(ee_candidate.edge0_index, 1)),
                        vertices.row(edges(ee_candidate.edge1_index, 0)),
                        vertices.row(edges(ee_candidate.edge1_index, 1)),
                        displacements.row(edges(ee_candidate.edge0_index, 0)),
                        displacements.row(edges(ee_candidate.edge0_index, 1)),
                        displacements.row(edges(ee_candidate.edge1_index, 0)),
                        displacements.row(edges(ee_candidate.edge1_index, 1)),
                        ee_candidate, ee_impacts);
                });
        },

        [&] {
            fv_impacts.clear();
            tbb::parallel_for_each(
                fv_candidates, [&](const FaceVertexCandidate& fv_candidate) {
                    detect_face_vertex_collisions_narrow_phase(
                        vertices.row(faces(fv_candidate.face_index, 0)),
                        vertices.row(faces(fv_candidate.face_index, 1)),
                        vertices.row(faces(fv_candidate.face_index, 2)),
                        vertices.row(fv_candidate.vertex_index),
                        displacements.row(faces(fv_candidate.face_index, 0)),
                        displacements.row(faces(fv_candidate.face_index, 1)),
                        displacements.row(faces(fv_candidate.face_index, 2)),
                        displacements.row(fv_candidate.vertex_index),
                        fv_candidate, fv_impacts);
                });
        });

    PROFILE_END(NARROW_PHASE);
    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
void detect_edge_vertex_collisions_narrow_phase(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const EdgeVertexCandidate& ev_candidate,
    EdgeVertexImpacts& ev_impacts)
{
    double toi, alpha;

    bool are_colliding = ccd::autodiff::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
    if (are_colliding) {
        ev_impacts.emplace_back(
            toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index);
    }
}

void detect_edge_edge_collisions_narrow_phase(
    const Eigen::VectorXd& Vi,
    const Eigen::VectorXd& Vj,
    const Eigen::VectorXd& Vk,
    const Eigen::VectorXd& Vl,
    const Eigen::VectorXd& Ui,
    const Eigen::VectorXd& Uj,
    const Eigen::VectorXd& Uk,
    const Eigen::VectorXd& Ul,
    const EdgeEdgeCandidate& ee_candidate,
    EdgeEdgeImpacts& ee_impacts)
{
    double toi;

    bool are_colliding = CTCD::edgeEdgeCTCD(
        Vi.head<3>(), Vj.head<3>(), Vk.head<3>(), Vl.head<3>(),
        Vi.head<3>() + Ui.head<3>(), Vj.head<3>() + Uj.head<3>(),
        Vk.head<3>() + Uk.head<3>(), Vl.head<3>() + Ul.head<3>(), /*eta=*/0,
        toi);
    if (are_colliding) {
        ee_impacts.emplace_back(
            toi, ee_candidate.edge0_index, -1, ee_candidate.edge1_index, -1);
    }
}

void detect_face_vertex_collisions_narrow_phase(
    const Eigen::VectorXd& Vi,
    const Eigen::VectorXd& Vj,
    const Eigen::VectorXd& Vk,
    const Eigen::VectorXd& Vl,
    const Eigen::VectorXd& Ui,
    const Eigen::VectorXd& Uj,
    const Eigen::VectorXd& Uk,
    const Eigen::VectorXd& Ul,
    const FaceVertexCandidate& fv_candidate,
    FaceVertexImpacts& fv_impacts)
{
    double toi;

    bool are_colliding = CTCD::vertexFaceCTCD(
        Vl.head<3>(), Vi.head<3>(), Vj.head<3>(), Vk.head<3>(),
        Vl.head<3>() + Ul.head<3>(), Vi.head<3>() + Ui.head<3>(),
        Vj.head<3>() + Uj.head<3>(), Vk.head<3>() + Uk.head<3>(), /*eta=*/0,
        toi);
    if (are_colliding) {
        fv_impacts.emplace_back(
            toi, fv_candidate.face_index, -1, fv_candidate.vertex_index, -1);
    }
}

} // namespace ccd
