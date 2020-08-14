// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "collision_detection.hpp"

#include <ccd/hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

void detect_collisions(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method)
{
    assert(vertices_t0.size() == vertices_t1.size());
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    // Do the broad phase by detecting candidate impacts
    Candidates candidates;
    detect_collision_candidates(
        vertices_t0, vertices_t1, edges, faces, group_ids, collision_types,
        candidates, method);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_collisions_from_candidates(
        vertices_t0, vertices_t1, edges, faces, candidates, impacts);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    assert(
        method == DetectionMethod::BRUTE_FORCE
        || method == DetectionMethod::HASH_GRID);
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    NAMED_PROFILE_POINT("collisions_detection__broad_phase", BROAD_PHASE);
    PROFILE_START(BROAD_PHASE);

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            vertices_t0, edges, faces, group_ids, collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_hash_grid(
            vertices_t0, vertices_t1, edges, faces, group_ids, collision_types,
            candidates, inflation_radius);
        break;
    }

    PROFILE_END(BROAD_PHASE);
}

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    std::vector<EdgeVertexCandidate>& ev_candidates)
{
    const bool check_group = group_ids.size() > 0;
    for (int ei = 0; ei < edges.rows(); ei++) {
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
}

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    std::vector<EdgeEdgeCandidate>& ee_candidates)
{
    const bool check_group = group_ids.size() > 0;
    for (int ei = 0; ei < edges.rows(); ei++) {
        // Loop over all remaining edges
        for (int ej = ei + 1; ej < edges.rows(); ej++) {
            bool has_common_endpoint = edges(ei, 0) == edges(ej, 0)
                || edges(ei, 0) == edges(ej, 1) || edges(ei, 1) == edges(ej, 0)
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

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    std::vector<FaceVertexCandidate>& fv_candidates)
{
    const bool check_group = group_ids.size() > 0;
    // Loop over all faces
    for (int fi = 0; fi < faces.rows(); fi++) {
        // Loop over all vertices
        for (int vi = 0; vi < vertices.rows(); vi++) {
            // Check that the vertex is not an endpoint of the edge
            bool is_endpoint =
                vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi, 2);
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

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all edges and all vertices.
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates)
{
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    // Loop over all edges
    tbb::parallel_invoke(
        [&] {
            if (collision_types & CollisionType::EDGE_VERTEX) {
                detect_edge_vertex_collision_candidates_brute_force(
                    vertices, edges, group_ids, candidates.ev_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::EDGE_EDGE) {
                detect_edge_edge_collision_candidates_brute_force(
                    edges, group_ids, candidates.ee_candidates);
            }
        },
        [&] {
            if (collision_types & CollisionType::FACE_VERTEX) {
                detect_face_vertex_collision_candidates_brute_force(
                    vertices, faces, group_ids, candidates.fv_candidates);
            }
        });
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_hash_grid(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius)
{
    using namespace CollisionType;
    HashGrid hashgrid;
    assert(edges.size()); // Even face-vertex need the edges
    hashgrid.resize(vertices_t0, vertices_t1, edges, inflation_radius);
    tbb::parallel_invoke(
        [&] {
            if (collision_types & (EDGE_VERTEX | FACE_VERTEX)) {
                hashgrid.addVertices(
                    vertices_t0, vertices_t1, inflation_radius);
            }
        },
        [&] {
            if (collision_types & (EDGE_VERTEX | EDGE_EDGE)) {
                hashgrid.addEdges(
                    vertices_t0, vertices_t1, edges, inflation_radius);
            }
        },
        [&] {
            if (collision_types & FACE_VERTEX) {
                hashgrid.addFaces(
                    vertices_t0, vertices_t1, faces, inflation_radius);
            }
        });

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    tbb::parallel_invoke(
        [&] {
            if (collision_types & EDGE_VERTEX) {
                hashgrid.getVertexEdgePairs(
                    edges, group_ids, candidates.ev_candidates);
            }
        },
        [&] {
            if (collision_types & EDGE_EDGE) {
                hashgrid.getEdgeEdgePairs(
                    edges, group_ids, candidates.ee_candidates);
            }
        },
        [&] {
            if (collision_types & FACE_VERTEX) {
                hashgrid.getFaceVertexPairs(
                    faces, group_ids, candidates.fv_candidates);
            }
        });
}

} // namespace ccd
