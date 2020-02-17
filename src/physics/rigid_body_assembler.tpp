#pragma once
#include "rigid_body_assembler.hpp"

namespace ccd {
namespace physics {

    template <typename T>
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    RigidBodyAssembler::world_vertices(const std::vector<Pose<T>>& poses) const
    {
        assert(poses.size() == num_bodies());
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> V(
            num_vertices(), dim());
        for (size_t i = 0; i < num_bodies(); ++i) {
            auto& rb = m_rbs[i];
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), rb.dim()) =
                rb.world_vertices(poses[i]);
        }
        return V;
    }

} // namespace physics
} // namespace ccd
