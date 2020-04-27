#pragma once
#include "distance_barrier_constraint.hpp"

#include <barrier/barrier.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const T& distance, const double eps) const
    {
        return barrier(distance - min_distance, eps, barrier_type);
    }

    template <typename T>
    void DistanceBarrierConstraint::compute_candidates_constraints(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<T>& poses,
        const Candidates& candidates,
        Eigen::VectorX<T>& barriers)
    {
        // distance barrier is evaluated at end-positions
        Eigen::MatrixX<T> vertices_t1 = bodies.world_vertices(poses);

        int num_barriers = candidates.size();
        barriers.resize(num_barriers);

        // Add edge-vertex barriers
        size_t start_i = 0;
        for (size_t i = 0; i < candidates.ev_candidates.size(); i++) {
            const auto& ev_candidate = candidates.ev_candidates[i];

            const Eigen::Vector2<T>& vertex =
                vertices_t1.row(ev_candidate.vertex_index);
            const Eigen::Vector2<T>& edge_v0 =
                vertices_t1.row(bodies.m_edges(ev_candidate.edge_index, 0));
            const Eigen::Vector2<T>& edge_v1 =
                vertices_t1.row(bodies.m_edges(ev_candidate.edge_index, 1));

            T distance = ccd::geometry::point_segment_distance<T>(
                vertex, edge_v0, edge_v1);

            barriers(start_i + i) = distance_barrier(distance);
        }

        // Add edge-edge barriers
        start_i += candidates.ev_candidates.size();
        for (size_t i = 0; i < candidates.ee_candidates.size(); i++) {
            const auto& ee_candidate = candidates.ee_candidates[i];

            const Eigen::Vector3<T>& edge0_v0 =
                vertices_t1.row(bodies.m_edges(ee_candidate.edge0_index, 0));
            const Eigen::Vector3<T>& edge0_v1 =
                vertices_t1.row(bodies.m_edges(ee_candidate.edge0_index, 1));
            const Eigen::Vector3<T>& edge1_v0 =
                vertices_t1.row(bodies.m_edges(ee_candidate.edge1_index, 0));
            const Eigen::Vector3<T>& edge1_v1 =
                vertices_t1.row(bodies.m_edges(ee_candidate.edge1_index, 1));

            T distance = ccd::geometry::segment_segment_distance<T>(
                edge0_v0, edge0_v1, edge1_v0, edge1_v1);

            barriers(start_i + i) = distance_barrier(distance);
        }

        // Add face-vertex barriers
        start_i += candidates.ee_candidates.size();
        for (size_t i = 0; i < candidates.fv_candidates.size(); i++) {
            const auto& fv_candidate = candidates.fv_candidates[i];

            const Eigen::Vector3<T>& vertex =
                vertices_t1.row(fv_candidate.vertex_index);
            const Eigen::Vector3<T>& face_v0 =
                vertices_t1.row(bodies.m_faces(fv_candidate.face_index, 0));
            const Eigen::Vector3<T>& face_v1 =
                vertices_t1.row(bodies.m_faces(fv_candidate.face_index, 1));
            const Eigen::Vector3<T>& face_v2 =
                vertices_t1.row(bodies.m_faces(fv_candidate.face_index, 2));

            T distance = ccd::geometry::point_triangle_distance<T>(
                vertex, face_v0, face_v1, face_v2);

            barriers(start_i + i) = distance_barrier(distance);
        }
    }
} // namespace opt
} // namespace ccd
