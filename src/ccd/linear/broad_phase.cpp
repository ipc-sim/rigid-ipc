// Detection collisions between different geometry.
// Includes continuous collision detection to compute the time of impact.
// Supported geometry: point vs edge

#include "broad_phase.hpp"

#include <tbb/parallel_invoke.h>

#include <ccd/rigid/broad_phase.hpp>
#include <logger.hpp>
#include <profiler.hpp>

using ipc::AABB;

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_linear(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
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
    ipc::Candidates& candidates,
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
    std::vector<ipc::EdgeVertexCandidate>& ev_candidates)
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
    std::vector<ipc::EdgeEdgeCandidate>& ee_candidates)
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
    std::vector<ipc::FaceVertexCandidate>& fv_candidates)
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
    ipc::Candidates& candidates)
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
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    using namespace CollisionType;
    ipc::HashGrid hashgrid;
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

    // Assume checking if vertex is and end-point of the edge is done by
    // `hashgrid.getVertexEdgePairs(...)`.
    if (collision_types & EDGE_VERTEX) {
        hashgrid.getVertexEdgePairs(edges, group_ids, candidates.ev_candidates);
    }
    if (collision_types & EDGE_EDGE) {
        hashgrid.getEdgeEdgePairs(edges, group_ids, candidates.ee_candidates);
    }
    if (collision_types & FACE_VERTEX) {
        hashgrid.getFaceVertexPairs(faces, group_ids, candidates.fv_candidates);
    }
}

void detect_collision_candidates_linear_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    if (bodies.dim() != 3 || (collision_types & CollisionType::EDGE_VERTEX)) {
        throw NotImplementedError(
            "detect_collision_candidates_rigid_bvh is not implemented for 2D "
            "or codimensional objects!");
    }

    std::vector<std::pair<int, int>> body_pairs =
        bodies.close_bodies(poses_t0, poses_t1, inflation_radius);

    typedef tbb::enumerable_thread_specific<ipc::Candidates> LocalStorage;
    LocalStorage storages;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), body_pairs.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            LocalStorage::reference loc_storage_candidates = storages.local();
            for (long i = range.begin(); i != range.end(); ++i) {
                int small_body_id = body_pairs[i].first;
                int large_body_id = body_pairs[i].second;
                if (bodies[small_body_id].num_faces()
                    > bodies[large_body_id].num_faces()) {
                    std::swap(small_body_id, large_body_id);
                }

                detect_collision_candidates_linear_bvh(
                    bodies, poses_t0, poses_t1, collision_types, small_body_id,
                    large_body_id, loc_storage_candidates, inflation_radius);
            }
        });

    merge_local_candidates(storages, candidates);
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_linear_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    const size_t bodyA_id,
    const size_t bodyB_id,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    const physics::RigidBody& bodyA = bodies[bodyA_id];
    const physics::RigidBody& bodyB = bodies[bodyB_id];

    // Compute the smaller body's vertices in the larger body's local
    // coordinates.
    const auto& pA_t0 = poses_t0[bodyA_id].position;
    const auto& pB_t0 = poses_t0[bodyB_id].position;
    const auto& pA_t1 = poses_t1[bodyA_id].position;
    const auto& pB_t1 = poses_t1[bodyB_id].position;
    const auto RA_t0 = poses_t0[bodyA_id].construct_rotation_matrix();
    const auto RB_t0 = poses_t0[bodyB_id].construct_rotation_matrix();
    const auto RA_t1 = poses_t1[bodyA_id].construct_rotation_matrix();
    const auto RB_t1 = poses_t1[bodyB_id].construct_rotation_matrix();

    // PROFILE_POINT("detect_collision_candidates_linear_bvh:compute_vertices");
    // PROFILE_START();
    const Eigen::MatrixXd VA_t0 =
        ((bodyA.vertices * RA_t0.transpose()).rowwise()
         + (pA_t0 - pB_t0).transpose())
        * RB_t0;
    const Eigen::MatrixXd VA_t1 =
        ((bodyA.vertices * RA_t1.transpose()).rowwise()
         + (pA_t1 - pB_t1).transpose())
        * RB_t1;
    // PROFILE_END();

    const Eigen::MatrixXd& VB = bodyB.vertices;
    const Eigen::MatrixXi &EA = bodyA.edges, &EB = bodyB.edges,
                          &FA = bodyA.faces, &FB = bodyB.faces;

    // For each face in the small body
    // std::mutex fv_mutex, ee_mutex;
    // tbb::parallel_for(0l, FA.rows(), [&](long fa_id) {
    for (long fa_id = 0; fa_id < FA.rows(); fa_id++) {
        // Construct a bbox of bodyA's face
        AABB fa_aabb = face_aabb(VA_t0, VA_t1, FA.row(fa_id), inflation_radius);

        std::vector<unsigned int> fb_ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            fa_aabb.getMin().array() - inflation_radius,
            fa_aabb.getMax().array() + inflation_radius, //
            fb_ids);

        for (const unsigned int fb_id : fb_ids) {
            AABB fb_aabb = face_aabb(VB, FB.row(fb_id), inflation_radius);

            for (int f_vi = 0; f_vi < 3; f_vi++) {
                // vertex A - face B
                long va_id = FA(fa_id, f_vi);
                if (bodyA.mesh_selector.vertex_to_face(va_id) == fa_id) {
                    AABB va_aabb =
                        vertex_aabb(VA_t0, VA_t1, va_id, inflation_radius);

                    if (AABB::are_overlaping(va_aabb, fb_aabb)) {
                        // std::scoped_lock lock(fv_mutex);
                        // Convert the local ids to the global ones
                        candidates.fv_candidates.emplace_back(
                            bodies.m_body_face_id[bodyB_id] + fb_id,
                            bodies.m_body_vertex_id[bodyA_id] + va_id);
                    }
                }

                // face A - vertex B
                long vb_id = FB(fb_id, f_vi);
                if (bodyB.mesh_selector.vertex_to_face(vb_id) == fb_id) {
                    AABB vb_aabb = vertex_aabb(VB, vb_id, inflation_radius);

                    if (AABB::are_overlaping(fa_aabb, vb_aabb)) {
                        // std::scoped_lock lock(fv_mutex);
                        // Convert the local ids to the global ones
                        candidates.fv_candidates.emplace_back(
                            bodies.m_body_face_id[bodyA_id] + fa_id,
                            bodies.m_body_vertex_id[bodyB_id] + vb_id);
                    }
                }
            }

            for (int fa_ei = 0; fa_ei < 3; fa_ei++) {
                long ea_id = bodyA.mesh_selector.face_to_edge(fa_id, fa_ei);

                if (bodyA.mesh_selector.edge_to_face(ea_id) != fa_id) {
                    continue;
                }

                AABB ea_aabb =
                    edge_aabb(VA_t0, VA_t1, EA.row(ea_id), inflation_radius);

                for (int fb_ei = 0; fb_ei < 3; fb_ei++) {
                    long eb_id = bodyB.mesh_selector.face_to_edge(fb_id, fb_ei);

                    if (bodyB.mesh_selector.edge_to_face(eb_id) != fb_id) {
                        continue;
                    }

                    // AABB edge-edge check
                    AABB eb_aabb =
                        edge_aabb(VB, EB.row(eb_id), inflation_radius);

                    if (AABB::are_overlaping(ea_aabb, eb_aabb)) {
                        // std::scoped_lock lock(ee_mutex);
                        // Convert the local ids to the global ones
                        candidates.ee_candidates.emplace_back(
                            bodies.m_body_edge_id[bodyA_id] + ea_id,
                            bodies.m_body_edge_id[bodyB_id] + eb_id);
                    }
                }
            }
        }
    }
    // );
}

} // namespace ccd
