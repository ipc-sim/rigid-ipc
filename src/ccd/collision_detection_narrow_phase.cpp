// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <autogen/time_of_impact_coeff.hpp>
#include <ccd/time_of_impact.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_edge_vertex_collisions_from_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const EdgeVertexCandidates& ev_candidates,
    EdgeVertexImpacts& ev_impacts)
{
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    ev_impacts.clear();

    for (const EdgeVertexCandidate& ev_candidate : ev_candidates) {
        detect_edge_vertex_collisions_narrow_phase(
            vertices.row(edges(ev_candidate.edge_index, 0)),
            vertices.row(edges(ev_candidate.edge_index, 1)),
            vertices.row(ev_candidate.vertex_index),
            displacements.row(edges(ev_candidate.edge_index, 0)),
            displacements.row(edges(ev_candidate.edge_index, 1)),
            displacements.row(ev_candidate.vertex_index), ev_candidate,
            ev_impacts);
    }

    PROFILE_END(NARROW_PHASE);
    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
void detect_edge_vertex_collisions_narrow_phase(const Eigen::Vector2d& Vi,
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
        ev_impacts.push_back(EdgeVertexImpact(
            toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index));
    }
}

// Compute the time of impact of a point and edge moVing in 2D.
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& /*Vi*/,
    const Eigen::Vector2d& /*Vj*/,
    const Eigen::Vector2d& /*Vk*/,
    const Eigen::Vector2d& /*Ui*/,
    const Eigen::Vector2d& /*Uj*/,
    const Eigen::Vector2d& /*Uk*/,
    double& /*toi*/,
    double& /*alpha*/,
    const double /*tolerance*/)
{
    throw DeprecatedError("use ccd::autodiff::compute_edge_vertex_time_of_impact instead");
}



// Convert a temporal parameterization to a spatial parameterization.
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& /*Vi*/,
    const Eigen::Vector2d& /*Vj*/,
    const Eigen::Vector2d& /*Vk*/,
    const Eigen::Vector2d& /*Ui*/,
    const Eigen::Vector2d& /*Uj*/,
    const Eigen::Vector2d& /*Uk*/,
    const double /*t*/,
    double& /*alpha*/,
    const double /*tolerance*/)
{
     throw DeprecatedError("use ccd::autodiff::temporal_parameterization_to_spatial instead");
}

} // namespace ccd
