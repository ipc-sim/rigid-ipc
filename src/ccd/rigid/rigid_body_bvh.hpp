#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/broad_phase/hash_grid.hpp>

#include <ccd/ccd.hpp>
#include <interval/interval.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

inline void
sort_body_pair(const RigidBodyAssembler& bodies, int& bodyA_id, int& bodyB_id)
{
    if (bodies[bodyA_id].bvh_size() > bodies[bodyB_id].bvh_size()) {
        std::swap(bodyA_id, bodyB_id);
    }
}

inline AABB vertex_aabb(const VectorMax3d& v, double inflation_radius = 0)
{
    return AABB(v.array() - inflation_radius, v.array() + inflation_radius);
}

inline AABB vertex_aabb(const VectorMax3I& v, double inflation_radius = 0)
{
    assert(v.rows() == 1 || v.cols() == 1);
    VectorMax3d min(v.size());
    VectorMax3d max(v.size());
    for (int i = 0; i < v.size(); i++) {
        if (empty(v(i))) {
            throw "interval is empty";
        }
        min(i) = (v(i) - inflation_radius).lower();
        max(i) = (v(i) + inflation_radius).upper();
    }
    return AABB(min, max);
}

template <typename T>
inline std::vector<AABB>
vertex_aabbs(const MatrixX<T>& V, double inflation_radius = 0)
{
    std::vector<AABB> aabbs;
    aabbs.reserve(V.rows());
    for (size_t i = 0; i < V.rows(); i++) {
        aabbs.push_back(vertex_aabb(VectorMax3<T>(V.row(i)), inflation_radius));
    }
    return aabbs;
}

void detect_body_pair_collision_candidates_from_aabbs(
    const RigidBodyAssembler& bodies,
    const std::vector<AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0);

template <typename T>
inline void detect_body_pair_collision_candidates_bvh(
    const RigidBodyAssembler& bodies,
    const Poses<T>& poses,
    int bodyA_id,
    int bodyB_id,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius = 0.0)
{
    sort_body_pair(bodies, bodyA_id, bodyB_id);

    // Compute the smaller body's vertices in the larger body's local
    // coordinates.
    // PROFILE_POINT(fmt::format(
    //     "detect_body_pair_collision_candidates_bvh<{}>:compute_vertices",
    //     get_type_name<T>()));
    // PROFILE_START();
    const auto RA = poses[bodyA_id].construct_rotation_matrix();
    const auto RB = poses[bodyB_id].construct_rotation_matrix();
    const auto& pA = poses[bodyA_id].position;
    const auto& pB = poses[bodyB_id].position;
    const MatrixX<T> VA =
        ((bodies[bodyA_id].vertices * RA.transpose()).rowwise()
         + (pA - pB).transpose())
        * RB;
    // PROFILE_END();

    detect_body_pair_collision_candidates_from_aabbs(
        bodies, vertex_aabbs(VA, inflation_radius), bodyA_id, bodyB_id,
        collision_types, candidates, inflation_radius);
}

void detect_body_pair_intersection_candidates_from_aabbs(
    const RigidBodyAssembler& bodies,
    const std::vector<AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    std::vector<EdgeFaceCandidate>& candidates,
    const double inflation_radius = 0.0);

} // namespace ipc::rigid
