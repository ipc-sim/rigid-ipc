#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>

#include <ccd/ccd.hpp>
#include <interval/interval.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

inline void sort_body_pair(
    const physics::RigidBodyAssembler& bodies, int& bodyA_id, int& bodyB_id)
{
    if (bodies[bodyA_id].bvh_size() > bodies[bodyB_id].bvh_size()) {
        std::swap(bodyA_id, bodyB_id);
    }
}

inline ipc::AABB
vertex_aabb(const Eigen::VectorX3d& v, double inflation_radius = 0)
{
    return ipc::AABB(
        v.array() - inflation_radius, v.array() + inflation_radius);
}

inline ipc::AABB
vertex_aabb(const Eigen::VectorX3I& v, double inflation_radius = 0)
{
    assert(v.rows() == 1 || v.cols() == 1);
    Eigen::VectorX3d min(v.size());
    Eigen::VectorX3d max(v.size());
    for (int i = 0; i < v.size(); i++) {
        if (empty(v(i))) {
            throw "interval is empty";
        }
        min(i) = (v(i) - inflation_radius).lower();
        max(i) = (v(i) + inflation_radius).upper();
    }
    return ipc::AABB(min, max);
}

template <typename T>
inline std::vector<ipc::AABB>
vertex_aabbs(const Eigen::MatrixX<T>& V, double inflation_radius = 0)
{
    std::vector<ipc::AABB> aabbs;
    aabbs.reserve(V.rows());
    for (size_t i = 0; i < V.rows(); i++) {
        aabbs.push_back(
            vertex_aabb(Eigen::VectorX3<T>(V.row(i)), inflation_radius));
    }
    return aabbs;
}

void detect_body_pair_collision_candidates_from_aabbs(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<ipc::AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    const int collision_types,
    ipc::Candidates& candidates,
    const double inflation_radius = 0.0);

template <typename T>
inline void detect_body_pair_collision_candidates_bvh(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<T>& poses,
    int bodyA_id,
    int bodyB_id,
    const int collision_types,
    ipc::Candidates& candidates,
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
    const Eigen::MatrixX<T> VA =
        ((bodies[bodyA_id].vertices * RA.transpose()).rowwise()
         + (pA - pB).transpose())
        * RB;
    // PROFILE_END();

    detect_body_pair_collision_candidates_from_aabbs(
        bodies, vertex_aabbs(VA, inflation_radius), bodyA_id, bodyB_id,
        collision_types, candidates, inflation_radius);
}

void detect_body_pair_intersection_candidates_from_aabbs(
    const physics::RigidBodyAssembler& bodies,
    const std::vector<ipc::AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    std::vector<ipc::EdgeFaceCandidate>& candidates,
    const double inflation_radius = 0.0);

} // namespace ccd
