// Detection collisions between different geometry.
// Includes continuous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "broad_phase.hpp"

#include <tbb/parallel_invoke.h>

#include <ccd/rigid/broad_phase.hpp>
#include <ccd/rigid/rigid_body_bvh.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace ipc::rigid {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_linear(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    PROFILE_POINT("detect_continous_collision_candidates_linear");
    PROFILE_START();

    switch (method) {
    case BRUTE_FORCE:
    case HASH_GRID: {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses_t0);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses_t1);
        detect_collision_candidates(
            V_t0, V_t1, bodies.m_edges, bodies.m_faces, bodies.group_ids(),
            collision_types, candidates, method, inflation_radius);
        break;
    }
    case BVH:
        detect_collision_candidates_linear_bvh(
            bodies, poses_t0, poses_t1, collision_types, candidates,
            inflation_radius);
        break;
    }

    PROFILE_END();
}

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
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    PROFILE_POINT("detect_collision_candidates(vertices...)");
    PROFILE_START();

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
    default:
        throw NotImplementedError(
            "detect_collision_candidates(vertices...) is only implemented for "
            "BRUTE_FORCE and HASH_GRID!");
    }

    PROFILE_END();
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

    if (collision_types & (EDGE_VERTEX | FACE_VERTEX)) {
        hashgrid.addVertices(vertices_t0, vertices_t1, inflation_radius);
    }

    if (collision_types & (EDGE_VERTEX | EDGE_EDGE)) {
        hashgrid.addEdges(vertices_t0, vertices_t1, edges, inflation_radius);
    }

    if (collision_types & FACE_VERTEX) {
        hashgrid.addFaces(vertices_t0, vertices_t1, faces, inflation_radius);
    }

    auto can_vertices_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids[vi] != group_ids[vj];
    };

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    if (collision_types & EDGE_VERTEX) {
        hashgrid.getVertexEdgePairs(
            edges, candidates.ev_candidates, can_vertices_collide);
    }
    if (collision_types & EDGE_EDGE) {
        hashgrid.getEdgeEdgePairs(
            edges, candidates.ee_candidates, can_vertices_collide);
    }
    if (collision_types & FACE_VERTEX) {
        hashgrid.getFaceVertexPairs(
            faces, candidates.fv_candidates, can_vertices_collide);
    }
}

void detect_collision_candidates_linear_bvh(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius)
{
    std::vector<std::pair<int, int>> body_pairs =
        bodies.close_bodies(poses_t0, poses_t1, inflation_radius);

    typedef tbb::enumerable_thread_specific<Candidates> LocalStorage;
    LocalStorage storages;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), body_pairs.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            LocalStorage::reference loc_storage_candidates = storages.local();
            for (long i = range.begin(); i != range.end(); ++i) {
                int bodyA_id = body_pairs[i].first;
                int bodyB_id = body_pairs[i].second;
                sort_body_pair(bodies, bodyA_id, bodyB_id);

                const auto& bodyA = bodies[bodyA_id];
                const auto& bodyB = bodies[bodyB_id];

                // PROFILE_POINT("detect_collision_candidates_linear_bvh:compute_vertices");
                // PROFILE_START();
                // Compute the smaller body's vertices in the larger body's
                // local coordinates.
                const auto& pA_t0 = poses_t0[bodyA_id].position;
                const auto& pB_t0 = poses_t0[bodyB_id].position;
                const auto& pA_t1 = poses_t1[bodyA_id].position;
                const auto& pB_t1 = poses_t1[bodyB_id].position;
                const auto RA_t0 =
                    poses_t0[bodyA_id].construct_rotation_matrix();
                const auto RB_t0 =
                    poses_t0[bodyB_id].construct_rotation_matrix();
                const auto RA_t1 =
                    poses_t1[bodyA_id].construct_rotation_matrix();
                const auto RB_t1 =
                    poses_t1[bodyB_id].construct_rotation_matrix();
                const Eigen::MatrixXd VA_t0 =
                    ((bodyA.vertices * RA_t0.transpose()).rowwise()
                     + (pA_t0 - pB_t0).transpose())
                    * RB_t0;
                const Eigen::MatrixXd VA_t1 =
                    ((bodyA.vertices * RA_t1.transpose()).rowwise()
                     + (pA_t1 - pB_t1).transpose())
                    * RB_t1;

                std::vector<AABB> VA_aabbs;
                VA_aabbs.reserve(bodyA.num_vertices());
                for (int i = 0; i < bodyA.num_vertices(); i++) {
                    const auto& v_t0 = VA_t0.row(i);
                    const auto& v_t1 = VA_t1.row(i);
                    VA_aabbs.emplace_back(
                        v_t0.cwiseMin(v_t1).array() - inflation_radius,
                        v_t0.cwiseMax(v_t1).array() + inflation_radius);
                }
                // PROFILE_END();

                detect_body_pair_collision_candidates_from_aabbs(
                    bodies, VA_aabbs, bodyA_id, bodyB_id, collision_types,
                    loc_storage_candidates, inflation_radius);
            }
        });

    merge_local_candidates(storages, candidates);
}

} // namespace ipc::rigid
