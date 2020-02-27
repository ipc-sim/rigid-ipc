#pragma once
#include "distance_barrier_constraint.hpp"

#include <barrier/barrier.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& a, // segment start
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b, // segment end
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c) // point
    {
        T distance = ccd::geometry::point_segment_distance<T>(c, a, b);
        return distance_barrier<T>(distance, m_barrier_epsilon);
    }

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const T distance, const double eps)
    {
        return barrier(distance - min_distance, eps, barrier_type);
    }

    template <typename T>
    void DistanceBarrierConstraint::compute_candidates_constraints(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        const physics::Poses<T>& displacements,
        const Candidates& candidates,
        Eigen::VectorX<T>& barriers)
    {
        // distance barrier is evaluated at end-positions
        physics::Poses<T> poses_t1 =
            physics::cast<double, T>(poses) + displacements;
        Eigen::MatrixX<T> vertices_t1 = bodies.world_vertices(poses_t1);

        int num_barriers = candidates.size();
        barriers.resize(num_barriers);
        barriers.setConstant(T(0.0));

        // Add edge-vertex barriers
        for (size_t i = 0; i < candidates.ev_candidates.size(); i++) {
            const auto& ev_candidate = candidates.ev_candidates[i];
            // a and b are the endpoints of the edge; c is the vertex
            long edge_id = ev_candidate.edge_index;
            long a_id = bodies.m_edges(edge_id, 0);
            long b_id = bodies.m_edges(edge_id, 1);
            long c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorX3<T> a = vertices_t1.row(a_id);
            Eigen::VectorX3<T> b = vertices_t1.row(b_id);
            Eigen::VectorX3<T> c = vertices_t1.row(c_id);

            barriers(int(i)) = distance_barrier<T>(a, b, c);
        }

        // TODO: Add edge-edge barriers
        // TODO: Add face-vertex barriers
        if (bodies.dim() != 2) {
            throw NotImplementedError(
                "DistanceBarrierConstraint::compute_candidates_constraints "
                "not implemented in 3D!");
        }
    }
} // namespace opt
} // namespace ccd
