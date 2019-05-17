// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include <ccd/collision_detection.hpp>
#include <ccd/hash.hpp>

#include <iostream>

#include <profiler.hpp>
#ifdef PROFILE_FUNCTIONS
long number_of_collision_detection_calls = 0;
double time_spent_detecting_collisions = 0;
#endif

#define EPSILON (1e-8)

namespace ccd {

// Find all edge-vertex collisions in one time step.
void detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    EdgeVertexImpacts& ev_impacts, DetectionMethod method, bool reset_impacts)
{
#ifdef PROFILE_FUNCTIONS
    number_of_collision_detection_calls++;
    igl::Timer timer;
    timer.start();
#endif

    if (reset_impacts) {
        ev_impacts.clear();
    }
    Eigen::MatrixXb skip_pair
        = Eigen::MatrixXb::Constant(vertices.rows(), edges.rows(), false);
    // If we do not recompute vertex impacts then we need to prevent duplicates
    for (EdgeVertexImpact ev_impact : ev_impacts) {
        skip_pair(ev_impact.vertex_index, ev_impact.edge_index) = true;
    }

    assert(vertices.size() == displacements.size());
    assert(method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_MAP);
    switch (method) {
    case BRUTE_FORCE:
        detect_edge_vertex_collisions_brute_force(
            vertices, displacements, edges, skip_pair, ev_impacts);
        break;
    case HASH_MAP:
        detect_edge_vertex_collisions_hash_map(
            vertices, displacements, edges, skip_pair, ev_impacts);
        break;
    }

#ifdef PROFILE_FUNCTIONS
    timer.stop();
    time_spent_detecting_collisions += timer.getElapsedTime();
#endif
}

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_edge_vertex_collisions_brute_force(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    const Eigen::MatrixXb& skip_pair, EdgeVertexImpacts& ev_impacts)
{
    // Loop over all edges
    for (int edge_id = 0; edge_id < edges.rows(); edge_id++) {
        // Loop over all vertices
        for (int vertex_id = 0; vertex_id < vertices.rows(); vertex_id++) {
            // Check that the vertex is not an endpoint of the edge
            if (!skip_pair(vertex_id, edge_id) && vertex_id != edges(edge_id, 0)
                && vertex_id != edges(edge_id, 1)) {
                // Check if the pair is colliding using the time of impact code.
                detect_edge_vertex_collisions_narrow_phase(
                    vertices.row(edges(edge_id, 0)),
                    vertices.row(edges(edge_id, 1)), vertices.row(vertex_id),
                    displacements.row(edges(edge_id, 0)),
                    displacements.row(edges(edge_id, 1)),
                    displacements.row(vertex_id), edge_id, vertex_id,
                    ev_impacts);
            }
        }
    }
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_edge_vertex_collisions_hash_map(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    const Eigen::MatrixXb& skip_pair, EdgeVertexImpacts& ev_impacts)
{
    Hash hashgrid;
    hashgrid.resize(vertices, displacements, edges);
    hashgrid.addVertices(vertices, displacements);
    hashgrid.addEdges(vertices, displacements, edges);

    Candidates candidates; // Set of vertex edge pairs
    hashgrid.getVertexEdgePairs(edges, candidates);

    std::cout << edges.rows() * vertices.rows() << " possible impacts\n"
              << candidates.size() << " candidate impacts" << std::endl;

    for (const Candidate& candidate : candidates) {
        const long vertex_id = candidate.first;
        const long edge_id = candidate.second;
        // Check if there is a collision between the vertex and edge. Assume
        // checking if vertex is and end-point of the edge is done by
        // `hashgrid.getVertexEdgePairs(...)`.
        if (!skip_pair(vertex_id, edge_id)) {
            // Check if the pair is colliding using the time of impact code
            detect_edge_vertex_collisions_narrow_phase(
                vertices.row(edges(edge_id, 0)),
                vertices.row(edges(edge_id, 1)), vertices.row(vertex_id),
                displacements.row(edges(edge_id, 0)),
                displacements.row(edges(edge_id, 1)),
                displacements.row(vertex_id), edge_id, vertex_id, ev_impacts);
        }
    }

    std::cout << ev_impacts.size() << " actual impacts" << std::endl;
}

} // namespace ccd
