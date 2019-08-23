#pragma once
#include "rigid_body_assembler.hpp"

namespace ccd {
namespace physics {

    template <typename T>
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    RigidBodyAssembler::world_vertices(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& positions) const
    {
        typedef Eigen::Matrix<T, 3, 1> Vector3d;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
        MatrixXd V(num_vertices(), 2);
        for (size_t i = 0; i < m_rbs.size(); ++i) {
            auto& rb = m_rbs[i];
            Vector3d p_i = positions.segment(3 * int(i), 3);
            V.block(m_body_vertex_id[i], 0, rb.vertices.rows(), 2)
                = rb.world_vertices(p_i);
        }
        return V;
    }

} // namespace physics
} // namespace ccd
