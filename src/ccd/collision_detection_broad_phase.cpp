// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <iostream>

#include <ccd/hash_grid.hpp>

#include <profiler.hpp>

namespace ccd {

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)

{
    return detect_edge_vertex_collisions(vertices, displacements, edges,
        Eigen::VectorXi(), ev_impacts, method);
}

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexImpacts& ev_impacts,
    DetectionMethod method)
{
    assert(vertices.size() == displacements.size());

    // Do the broad phase by detecting candidate impacts
    EdgeVertexCandidates ev_candidates;
    detect_edge_vertex_collision_candidates(
        vertices, displacements, edges, group_ids, ev_candidates, method);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_edge_vertex_collisions_from_candidates(vertices, displacements,
        edges, ev_candidates, ev_impacts);
}

void detect_edge_vertex_collision_candidates(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    assert(method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);

    PROFILE_START();
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_edge_vertex_collision_candidates_brute_force(
            vertices, displacements, edges, group_ids, ev_candidates);
        break;
    case HASH_GRID:
        detect_edge_vertex_collision_candidates_hash_grid(
            vertices, displacements, edges, group_ids, ev_candidates, inflation_radius);
        break;
    }

    PROFILE_END(BROAD_PHASE);
    PROFILE_END();
}

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& /* displacements */,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    const double /* inflation_radius */)
{
    const bool check_group = group_ids.size() > 0;
    // Loop over all edges
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        // Loop over all vertices
        for (int vertex_index = 0; vertex_index < vertices.rows();
             vertex_index++) {
            // Check that the vertex is not an endpoint of the edge
            bool is_endpoint = vertex_index == edges(edge_index, 0)
                || vertex_index == edges(edge_index, 1);
            bool same_group = false;
            if (check_group) {
                // TODO: Check for the other vertex of the edge too.
                same_group = group_ids(vertex_index)
                    == group_ids(edges(edge_index, 0));
            }
            if (!is_endpoint && !same_group) {
                ev_candidates.push_back(
                    EdgeVertexCandidate(edge_index, vertex_index));
            }
        }
    }
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_edge_vertex_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates,
    const double inflation_radius)
{
    HashGrid hashgrid;
    hashgrid.resize(vertices, displacements, edges, inflation_radius);
    hashgrid.addVertices(vertices, displacements, inflation_radius);
    hashgrid.addEdges(vertices, displacements, edges, inflation_radius);

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    hashgrid.getVertexEdgePairs(edges, group_ids, ev_candidates);
}

} // namespace ccd
