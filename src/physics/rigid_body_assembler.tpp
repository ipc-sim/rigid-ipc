#pragma once
#include "rigid_body_assembler.hpp"

namespace ccd {
namespace physics {

    template <typename T>
    Eigen::MatrixX<T>
    RigidBodyAssembler::world_vertices(const Poses<T>& poses) const
    {
        assert(poses.size() == num_bodies());
        Eigen::MatrixX<T> V(num_vertices(), dim());
        for (size_t i = 0; i < num_bodies(); ++i) {
            const RigidBody& rb = m_rbs[i];
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), rb.dim()) =
                rb.world_vertices(poses[i]);
        }
        return V;
    }

} // namespace physics
} // namespace ccd
