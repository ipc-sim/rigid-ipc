// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <ccd/hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

void detect_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexImpacts& ev_impacts,
    EdgeEdgeImpacts& ee_impacts,
    FaceVertexImpacts& fv_impacts,
    DetectionMethod method)
{
    assert(vertices.size() == displacements.size());

    // Do the broad phase by detecting candidate impacts
    EdgeVertexCandidates ev_candidates;
    EdgeEdgeCandidates ee_candidates;
    FaceVertexCandidates fv_candidates;
    detect_collision_candidates(
        vertices, displacements, edges, faces, group_ids, collision_types,
        ev_candidates, ee_candidates, fv_candidates, method);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_collisions_from_candidates(
        vertices, displacements, edges, faces, ev_candidates, ee_candidates,
        fv_candidates, ev_impacts, ee_impacts, fv_impacts);
}

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)
{
    // Both of these will be untouched
    EdgeEdgeImpacts ee_impacts;
    FaceVertexImpacts fv_impacts;
    detect_collisions(
        vertices, displacements, edges, Eigen::MatrixXi(), group_ids,
        CollisionType::EDGE_VERTEX, ev_impacts, ee_impacts, fv_impacts, method);
}

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)
{
    return detect_edge_vertex_collisions(
        vertices, displacements, edges, Eigen::VectorXi(), ev_impacts, method);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    assert(
        method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);

    PROFILE_START();
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            vertices, edges, faces, group_ids, collision_types, ev_candidates,
            ee_candidates, fv_candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_hash_grid(
            vertices, displacements, edges, faces, group_ids, collision_types,
            ev_candidates, ee_candidates, fv_candidates, inflation_radius);
        spdlog::debug(
            "hash_grid_ev_candidates.size()={:d} "
            "hash_grid_ee_candidates.size()={:d} "
            "hash_grid_fv_candidates.size()={:d}",
            ev_candidates.size(), ee_candidates.size(), fv_candidates.size());
        break;
    }

    PROFILE_END(BROAD_PHASE);
    PROFILE_END();
}

void detect_edge_vertex_collision_candidates(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    // Both of these will be untouched
    EdgeEdgeCandidates ee_candidates;
    FaceVertexCandidates fv_candidates;
    detect_collision_candidates(
        vertices, displacements, edges, Eigen::MatrixXi(), group_ids,
        CollisionType::EDGE_VERTEX, ev_candidates, ee_candidates, fv_candidates,
        method, inflation_radius);
}

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    const bool check_group = group_ids.size() > 0;
    // Loop over all edges
    for (int ei = 0; ei < edges.rows(); ei++) {
        if (collision_types & CollisionType::EDGE_VERTEX) {
            // Loop over all vertices
            for (int vi = 0; vi < vertices.rows(); vi++) {
                // Check that the vertex is not an endpoint of the edge
                bool is_endpoint = vi == edges(ei, 0) || vi == edges(ei, 1);
                bool same_group = check_group
                    && (group_ids(vi) == group_ids(edges(ei, 0))
                        || group_ids(vi) == group_ids(edges(ei, 1)));
                if (!is_endpoint && !same_group) {
                    ev_candidates.emplace_back(ei, vi);
                }
            }
        }
        if (collision_types & CollisionType::EDGE_EDGE) {
            // Loop over all remaining edges
            for (int ej = ei + 1; ej < edges.rows(); ej++) {
                bool has_common_endpoint = edges(ei, 0) == edges(ej, 0)
                    || edges(ei, 0) == edges(ej, 1)
                    || edges(ei, 1) == edges(ej, 0)
                    || edges(ei, 1) == edges(ej, 1);
                bool same_group = check_group
                    && (group_ids(edges(ei, 0)) == group_ids(edges(ej, 0))
                        || group_ids(edges(ei, 0)) == group_ids(edges(ej, 1))
                        || group_ids(edges(ei, 1)) == group_ids(edges(ej, 0))
                        || group_ids(edges(ei, 1)) == group_ids(edges(ej, 1)));
                if (!has_common_endpoint && !same_group) {
                    ee_candidates.emplace_back(ei, ej);
                }
            }
        }
    }
    if (collision_types & CollisionType::FACE_VERTEX) {
        // Loop over all faces
        for (int fi = 0; fi < faces.rows(); fi++) {
            // Loop over all vertices
            for (int vi = 0; vi < vertices.rows(); vi++) {
                // Check that the vertex is not an endpoint of the edge
                bool is_endpoint = vi == faces(fi, 0) || vi == faces(fi, 1)
                    || vi == faces(fi, 2);
                bool same_group = check_group
                    && (group_ids(vi) == group_ids(faces(fi, 0))
                        || group_ids(vi) == group_ids(faces(fi, 1))
                        || group_ids(vi) == group_ids(faces(fi, 2)));
                if (!is_endpoint && !same_group) {
                    fv_candidates.emplace_back(fi, vi);
                }
            }
        }
    }
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    EdgeVertexCandidates& ev_candidates,
    EdgeEdgeCandidates& ee_candidates,
    FaceVertexCandidates& fv_candidates,
    const double inflation_radius)
{
    using namespace CollisionType;
    HashGrid hashgrid;
    assert(edges.size()); // Even face-vertex need the edges
    hashgrid.resize(vertices, displacements, edges, inflation_radius);
    tbb::parallel_invoke(
        [&] {
            if (collision_types & (EDGE_VERTEX | FACE_VERTEX)) {
                hashgrid.addVertices(vertices, displacements, inflation_radius);
            }
        },
        [&] {
            if (collision_types & (EDGE_VERTEX | EDGE_EDGE)) {
                hashgrid.addEdges(
                    vertices, displacements, edges, inflation_radius);
            }
        },
        [&] {
            if (collision_types & FACE_VERTEX) {
                hashgrid.addFaces(
                    vertices, displacements, faces, inflation_radius);
            }
        });

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    tbb::parallel_invoke(
        [&] {
            if (collision_types & EDGE_VERTEX) {
                hashgrid.getVertexEdgePairs(edges, group_ids, ev_candidates);
            }
        },
        [&] {
            if (collision_types & EDGE_EDGE) {
                hashgrid.getEdgeEdgePairs(edges, group_ids, ee_candidates);
            }
        },
        [&] {
            if (collision_types & FACE_VERTEX) {
                hashgrid.getFaceVertexPairs(faces, group_ids, fv_candidates);
            }
        });
}

} // namespace ccd
