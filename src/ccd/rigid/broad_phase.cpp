#include "broad_phase.hpp"

#include <tbb/parallel_for.h>

#include <ccd/linear/broad_phase.hpp>
#include <ccd/rigid/rigid_body_hash_grid.hpp>
#include <logger.hpp>
#include <profiler.hpp>
#include <utils/type_name.hpp>

using ipc::AABB;

namespace ccd {

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Discrete Collision Detection
// NOTE: Yes, this is inside the CCD directory.
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_rigid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method,
    const double inflation_radius)
{
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    PROFILE_POINT("detect_discrete_collision_candidates_rigid");
    PROFILE_START();

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            bodies.world_vertices(poses), bodies.m_edges, bodies.m_faces,
            bodies.group_ids(), collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_rigid_hash_grid(
            bodies, poses, collision_types, candidates, inflation_radius);
        break;
    case BVH:
        detect_collision_candidates_rigid_bvh(
            bodies, poses, collision_types, candidates, inflation_radius);
        break;
    }

    PROFILE_END();
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    std::vector<std::pair<int, int>> body_pairs =
        bodies.close_bodies(poses, poses, inflation_radius);

    if (body_pairs.size() == 0) {
        return;
    }

    RigidBodyHashGrid hashgrid;
    hashgrid.resize(bodies, poses, body_pairs, inflation_radius);
    hashgrid.addBodies(bodies, poses, body_pairs, inflation_radius);

    const Eigen::VectorXi& group_ids = bodies.group_ids();
    if (collision_types & CollisionType::EDGE_VERTEX) {
        hashgrid.getVertexEdgePairs(
            bodies.m_edges, group_ids, candidates.ev_candidates);
    }
    if (collision_types & CollisionType::EDGE_EDGE) {
        hashgrid.getEdgeEdgePairs(
            bodies.m_edges, group_ids, candidates.ee_candidates);
    }
    if (collision_types & CollisionType::FACE_VERTEX) {
        hashgrid.getFaceVertexPairs(
            bodies.m_faces, group_ids, candidates.fv_candidates);
    }
}

// Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
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
        bodies.close_bodies(poses, poses, inflation_radius);

    // auto posesI = physics::cast<Interval>(poses);

    // for (const auto& body_pair : body_pairs) {
    ThreadSpecificCandidates storages;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), body_pairs.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            ThreadSpecificCandidates::reference local_storage_candidates =
                storages.local();
            for (long i = range.begin(); i != range.end(); ++i) {
                int small_body_id = body_pairs[i].first;
                int large_body_id = body_pairs[i].second;
                if (bodies[small_body_id].num_faces()
                    > bodies[large_body_id].num_faces()) {
                    std::swap(small_body_id, large_body_id);
                }

                detect_collision_candidates_rigid_bvh(
                    // bodies, posesI, collision_types, small_body_id,
                    bodies, poses, collision_types, small_body_id,
                    large_body_id, local_storage_candidates, inflation_radius);
            }
        });
    //}

    merge_local_candidates(storages, candidates);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Continous Collision Detection
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates_rigid(
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

    PROFILE_POINT("detect_continuous_collision_candidates_rigid");
    PROFILE_START();

    switch (method) {
    case BRUTE_FORCE:
        detect_collision_candidates_brute_force(
            bodies.world_vertices(poses_t0), bodies.m_edges, bodies.m_faces,
            bodies.group_ids(), collision_types, candidates);
        break;
    case HASH_GRID:
        detect_collision_candidates_rigid_hash_grid(
            bodies, poses_t0, poses_t1, collision_types, candidates,
            inflation_radius);
        break;
    case BVH:
        detect_collision_candidates_rigid_bvh(
            bodies, poses_t0, poses_t1, collision_types, candidates,
            inflation_radius);
        break;
    }

    PROFILE_END();
}

// Find all edge-vertex collisions in one time step using spatial-hashing to
// only compare points and edge in the same cells.
void detect_collision_candidates_rigid_hash_grid(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    std::vector<std::pair<int, int>> body_pairs =
        bodies.close_bodies(poses_t0, poses_t1, inflation_radius);

    if (body_pairs.size() == 0) {
        return;
    }

    RigidBodyHashGrid hashgrid;
    hashgrid.resize(bodies, poses_t0, poses_t1, body_pairs, inflation_radius);
    hashgrid.addBodies(
        bodies, poses_t0, poses_t1, body_pairs, inflation_radius);

    const Eigen::VectorXi& group_ids = bodies.group_ids();
    if (collision_types & CollisionType::EDGE_VERTEX) {
        hashgrid.getVertexEdgePairs(
            bodies.m_edges, group_ids, candidates.ev_candidates);
    }
    if (collision_types & CollisionType::EDGE_EDGE) {
        hashgrid.getEdgeEdgePairs(
            bodies.m_edges, group_ids, candidates.ee_candidates);
    }
    if (collision_types & CollisionType::FACE_VERTEX) {
        hashgrid.getFaceVertexPairs(
            bodies.m_faces, group_ids, candidates.fv_candidates);
    }
}

template <typename Derived>
inline AABB
intervals_to_AABB(const Eigen::MatrixBase<Derived>& x, double inflation_radius)
{
    assert(x.rows() == 1 || x.cols() == 1);
    Eigen::VectorX3d min(x.size());
    Eigen::VectorX3d max(x.size());
    for (int i = 0; i < x.size(); i++) {
        if (empty(x(i))) {
            throw "interval is empty";
        }
        min(i) = x(i).lower() - inflation_radius;
        max(i) = x(i).upper() + inflation_radius;
    }
    return AABB(min, max);
}

// Use a BVH to create a set of all candidate collisions.
void detect_collision_candidates_rigid_bvh(
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

    physics::Poses<Interval> poses = physics::interpolate(
        physics::cast<Interval>(poses_t0), physics::cast<Interval>(poses_t1),
        Interval(0, 1));

    ThreadSpecificCandidates storages;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), body_pairs.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            ThreadSpecificCandidates::reference local_storage_candidates =
                storages.local();
            for (long i = range.begin(); i != range.end(); ++i) {
                int small_body_id = body_pairs[i].first;
                int large_body_id = body_pairs[i].second;
                if (bodies[small_body_id].num_faces()
                    > bodies[large_body_id].num_faces()) {
                    std::swap(small_body_id, large_body_id);
                }

                detect_collision_candidates_rigid_bvh(
                    bodies, poses, collision_types, small_body_id,
                    large_body_id, local_storage_candidates, inflation_radius);
            }
        });

    merge_local_candidates(storages, candidates);
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

template <typename T>
void detect_collision_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<T>& poses,
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
    const auto RA = poses[bodyA_id].construct_rotation_matrix();
    const auto RB = poses[bodyB_id].construct_rotation_matrix();
    const auto& pA = poses[bodyA_id].position;
    const auto& pB = poses[bodyB_id].position;

    // PROFILE_POINT(fmt::format(
    //     "detect_collision_candidates_rigid_bvh<{}>:compute_vertices",
    //     get_type_name<T>()));
    // PROFILE_START();
    const Eigen::MatrixX<T> VA =
        ((bodyA.vertices * RA.transpose()).rowwise() + (pA - pB).transpose())
        * RB;
    // PROFILE_END();

    const Eigen::MatrixXd& VB = bodyB.vertices;
    const Eigen::MatrixXi &EA = bodyA.edges, &EB = bodyB.edges,
                          &FA = bodyA.faces, &FB = bodyB.faces;

    // For each face in the small body
    // std::mutex fv_mutex, ee_mutex;
    // tbb::parallel_for(0l, FA.rows(), [&](long fa_id) {
    for (long fa_id = 0; fa_id < FA.rows(); fa_id++) {
        // Construct a bbox of bodyA's face
        AABB fa_aabb = face_aabb(VA, FA.row(fa_id), inflation_radius);

        // DEBUG: print fa_aabb size vs size of triangle
        // if (std::is_same<T, Interval>::value) {
        //     double fa_longest_edge_length = (bodyA.vertices.row(FA(fa_id, 0))
        //                                      - bodyA.vertices.row(FA(fa_id,
        //                                      1)))
        //                                         .norm();
        //     fa_longest_edge_length = std::max(
        //         fa_longest_edge_length,
        //         (bodyA.vertices.row(FA(fa_id, 1))
        //          - bodyA.vertices.row(FA(fa_id, 2)))
        //             .norm());
        //     fa_longest_edge_length = std::max(
        //         fa_longest_edge_length,
        //         (bodyA.vertices.row(FA(fa_id, 2))
        //          - bodyA.vertices.row(FA(fa_id, 0)))
        //             .norm());
        //
        //     fmt::print(
        //         "face_aabb_diag_norm,{:g},fa_longest_edge_length,{:g}\n",
        //         (fa_aabb.getMax() - fa_aabb.getMin()).norm(),
        //         fa_longest_edge_length);
        // }

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
                    AABB va_aabb = vertex_aabb(VA, va_id, inflation_radius);

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

                AABB ea_aabb = edge_aabb(VA, EA.row(ea_id), inflation_radius);

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

void merge_local_candidates(
    const ThreadSpecificCandidates& storages, ipc::Candidates& candidates)
{
    PROFILE_POINT("merge_local_candidates");
    PROFILE_START();
    // size up the candidates
    size_t ev_size = 0, ee_size = 0, fv_size = 0;
    for (const auto& local_candidates : storages) {
        ev_size += local_candidates.ev_candidates.size();
        ee_size += local_candidates.ee_candidates.size();
        fv_size += local_candidates.fv_candidates.size();
    }
    // serial merge!
    candidates.ev_candidates.reserve(ev_size);
    candidates.ee_candidates.reserve(ee_size);
    candidates.fv_candidates.reserve(fv_size);
    for (const auto& local_candidates : storages) {
        candidates.ev_candidates.insert(
            candidates.ev_candidates.end(),
            local_candidates.ev_candidates.begin(),
            local_candidates.ev_candidates.end());
        candidates.ee_candidates.insert(
            candidates.ee_candidates.end(),
            local_candidates.ee_candidates.begin(),
            local_candidates.ee_candidates.end());
        candidates.fv_candidates.insert(
            candidates.fv_candidates.end(),
            local_candidates.fv_candidates.begin(),
            local_candidates.fv_candidates.end());
    }
    PROFILE_END();
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase Intersection Detection
///////////////////////////////////////////////////////////////////////////////

template <typename T>
void detect_intersection_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<T>& poses,
    const size_t bodyA_id,
    const size_t bodyB_id,
    std::vector<ipc::EdgeFaceCandidate>& ef_candidates)
{
    const physics::RigidBody& bodyA = bodies[bodyA_id];
    const physics::RigidBody& bodyB = bodies[bodyB_id];

    // Compute the smaller body's vertices in the larger body's local
    // coordinates.
    const auto RA = poses[bodyA_id].construct_rotation_matrix();
    const auto RB = poses[bodyB_id].construct_rotation_matrix();
    const auto& pA = poses[bodyA_id].position;
    const auto& pB = poses[bodyB_id].position;

    // PROFILE_POINT(fmt::format(
    //     "detect_intersection_candidates_rigid_bvh<{}>:compute_vertices",
    //     get_type_name<T>()));
    // PROFILE_START();
    const Eigen::MatrixX<T> VA =
        ((bodyA.vertices * RA.transpose()).rowwise() + (pA - pB).transpose())
        * RB;
    // PROFILE_END();

    const Eigen::MatrixXd& VB = bodyB.vertices;
    const Eigen::MatrixXi &EA = bodyA.edges, &EB = bodyB.edges,
                          &FA = bodyA.faces, &FB = bodyB.faces;

    // For each face in the small body
    // std::mutex fv_mutex, ee_mutex;
    // tbb::parallel_for(0l, FA.rows(), [&](long fa_id) {
    for (long fa_id = 0; fa_id < FA.rows(); fa_id++) {
        // Construct a bbox of bodyA's face
        AABB fa_aabb = face_aabb(VA, FA.row(fa_id));

        std::vector<unsigned int> fb_ids;
        bodyB.bvh.intersect_box(
            fa_aabb.getMin().array(), fa_aabb.getMax().array(), fb_ids);

        for (const unsigned int fb_id : fb_ids) {
            AABB fb_aabb = face_aabb(VB, FB.row(fb_id));

            for (int fa_ei = 0; fa_ei < 3; fa_ei++) {
                long ea_id = bodyA.mesh_selector.face_to_edge(fa_id, fa_ei);

                if (bodyA.mesh_selector.edge_to_face(ea_id) != fa_id) {
                    continue;
                }

                AABB ea_aabb = edge_aabb(VA, EA.row(ea_id));
                if (AABB::are_overlaping(ea_aabb, fb_aabb)) {
                    // std::scoped_lock lock(ee_mutex);
                    // Convert the local ids to the global ones
                    ef_candidates.emplace_back(
                        bodies.m_body_edge_id[bodyA_id] + ea_id,
                        bodies.m_body_face_id[bodyB_id] + fb_id);
                }
            }

            for (int fb_ei = 0; fb_ei < 3; fb_ei++) {
                long eb_id = bodyB.mesh_selector.face_to_edge(fb_id, fb_ei);

                if (bodyB.mesh_selector.edge_to_face(eb_id) != fb_id) {
                    continue;
                }

                AABB eb_aabb = edge_aabb(VB, EB.row(eb_id));
                if (AABB::are_overlaping(fa_aabb, eb_aabb)) {
                    // std::scoped_lock lock(ee_mutex);
                    // Convert the local ids to the global ones
                    ef_candidates.emplace_back(
                        bodies.m_body_edge_id[bodyB_id] + eb_id,
                        bodies.m_body_face_id[bodyA_id] + fa_id);
                }
            }
        }
    }
    // );
}

typedef tbb::enumerable_thread_specific<std::vector<ipc::EdgeFaceCandidate>>
    ThreadSpecificEFCandidates;

void merge_local_candidates(
    const ThreadSpecificEFCandidates& storages,
    std::vector<ipc::EdgeFaceCandidate>& ef_candidates)
{
    PROFILE_POINT("merge_local_ef_candidates");
    PROFILE_START();
    // size up the candidates
    size_t size = 0;
    for (const auto& local_candidates : storages) {
        size += local_candidates.size();
    }
    // serial merge
    ef_candidates.reserve(size);
    for (const auto& local_candidates : storages) {
        ef_candidates.insert(
            ef_candidates.end(), //
            local_candidates.begin(), local_candidates.end());
    }
    PROFILE_END();
}

// Use a BVH to create a set of all candidate intersections.
void detect_intersection_candidates_rigid_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    std::vector<ipc::EdgeFaceCandidate>& ef_candidates)
{
    std::vector<std::pair<int, int>> body_pairs =
        bodies.close_bodies(poses, poses, /*inflation_radius=*/0);

    auto posesI = physics::cast<Interval>(poses);

    // for (const auto& body_pair : body_pairs) {
    ThreadSpecificEFCandidates storages;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), body_pairs.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            ThreadSpecificEFCandidates::reference local_storage_candidates =
                storages.local();
            for (long i = range.begin(); i != range.end(); ++i) {
                int small_body_id = body_pairs[i].first;
                int large_body_id = body_pairs[i].second;
                if (bodies[small_body_id].num_faces()
                    > bodies[large_body_id].num_faces()) {
                    std::swap(small_body_id, large_body_id);
                }

                detect_intersection_candidates_rigid_bvh(
                    bodies, posesI, small_body_id, large_body_id,
                    local_storage_candidates);
            }
        });
    //}

    merge_local_candidates(storages, ef_candidates);
}

} // namespace ccd
