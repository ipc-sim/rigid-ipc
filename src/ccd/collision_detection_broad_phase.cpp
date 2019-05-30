// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <ccd/collision_detection.hpp>
#include <ccd/hash_grid.hpp>

#include <iostream>

#include <profiler.hpp>

namespace ccd {

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts, DetectionMethod method, bool reset_impacts)
{
    assert(vertices.size() == displacements.size());
    assert(method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);

    Eigen::MatrixXb skip_pair;
    EdgeVertexCandidates ev_candidates;
    PROFILE(
        if (reset_impacts) { ev_impacts.clear(); }

        skip_pair
        = Eigen::MatrixXb::Zero(edges.rows(), vertices.rows());
        // If we do not reset impacts then we need to prevent duplicates
        for (EdgeVertexImpact ev_impact
             : ev_impacts) {
            skip_pair(ev_impact.edge_index, ev_impact.vertex_index) = true;
        }

        switch (method) {
            case BRUTE_FORCE:
                detect_edge_vertex_collisions_brute_force(
                    vertices, displacements, edges, ev_candidates);
                break;
            case HASH_GRID:
                detect_edge_vertex_collisions_hash_map(
                    vertices, displacements, edges, ev_candidates);
                break;
        },
        ProfiledPoint::DETECTING_COLLISIONS_BROAD_PHASE);

    for (const EdgeVertexCandidate& ev_candidate : ev_candidates) {
        if (!skip_pair(ev_candidate.edge_index, ev_candidate.vertex_index)) {
            // Check if the pair is colliding using the time of impact code
            detect_edge_vertex_collisions_narrow_phase(
                vertices.row(edges(ev_candidate.edge_index, 0)),
                vertices.row(edges(ev_candidate.edge_index, 1)),
                vertices.row(ev_candidate.vertex_index),
                displacements.row(edges(ev_candidate.edge_index, 0)),
                displacements.row(edges(ev_candidate.edge_index, 1)),
                displacements.row(ev_candidate.vertex_index), ev_candidate,
                ev_impacts);
        }
    }
}

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_edge_vertex_collisions_brute_force(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexCandidates& ev_candidates)
{
    // Loop over all edges
    for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
        // Loop over all vertices
        for (int vertex_index = 0; vertex_index < vertices.rows();
             vertex_index++) {
            // Check that the vertex is not an endpoint of the edge
            if (vertex_index != edges(edge_index, 0)
                && vertex_index != edges(edge_index, 1)) {
                ev_candidates.push_back(
                    EdgeVertexCandidate(edge_index, vertex_index));
            }
        }
    }
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_edge_vertex_collisions_hash_map(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexCandidates& ev_candidates)
{
    HashGrid hashgrid;
    hashgrid.resize(vertices, displacements, edges);
    hashgrid.addVertices(vertices, displacements);
    hashgrid.addEdges(vertices, displacements, edges);

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    hashgrid.getVertexEdgePairs(edges, ev_candidates);
}

} // namespace ccd
