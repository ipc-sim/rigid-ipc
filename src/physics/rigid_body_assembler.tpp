#pragma once
#include "rigid_body_assembler.hpp"

namespace ipc::rigid {

template <typename T>
MatrixX<T> RigidBodyAssembler::world_vertices(const Poses<T>& poses) const
{
    assert(poses.size() == num_bodies());
    MatrixX<T> V(num_vertices(), dim());
    for (size_t i = 0; i < num_bodies(); ++i) {
        const RigidBody& rb = m_rbs[i];
        V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), rb.dim()) =
            rb.world_vertices(poses[i]);
    }
    return V;
}

template <typename T>
MatrixX<T> RigidBodyAssembler::world_vertices(
    const std::vector<MatrixMax3<T>>& rotations,
    const std::vector<VectorMax3<T>>& positions) const
{
    assert(rotations.size() == num_bodies());
    assert(positions.size() == num_bodies());
    MatrixX<T> V(num_vertices(), dim());
    for (size_t i = 0; i < num_bodies(); ++i) {
        const RigidBody& rb = m_rbs[i];
        V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), rb.dim()) =
            rb.world_vertices(rotations[i], positions[i]);
    }
    return V;
}

template <typename T>
VectorMax3<T> RigidBodyAssembler::world_vertex(
    const Pose<T>& pose, const int vertex_idx) const
{
    long body_idx, vertex_local_idx;
    global_to_local_vertex(vertex_idx, body_idx, vertex_local_idx);
    return m_rbs[body_idx].world_vertex(pose, vertex_local_idx);
}

template <typename T>
VectorMax3<T> RigidBodyAssembler::world_vertex(
    const Poses<T>& poses, const int vertex_idx) const
{
    assert(poses.size() == num_bodies());
    long body_idx, vertex_local_idx;
    global_to_local_vertex(vertex_idx, body_idx, vertex_local_idx);
    return m_rbs[body_idx].world_vertex(poses[body_idx], vertex_local_idx);
}

} // namespace ipc::rigid
