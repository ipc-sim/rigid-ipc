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
        if (use_log_barrier) {
            return opt::poly_log_barrier<T>(distance - min_distance, eps);
        } else {
            return opt::spline_barrier<T>(distance - min_distance, eps);
        }
    }

    template <typename T>
    void DistanceBarrierConstraint::compute_candidates_constraints(
        const Eigen::MatrixX<T>& Uk,
        const EdgeVertexCandidates& ev_candidates,
        Eigen::VectorX<T>& barriers)
    {
        // distance barrier is evaluated at end-positions
        Eigen::MatrixX<T> vertices_t1 = vertices.cast<T>() + Uk;

        barriers.resize(ev_candidates.size(), 1);
        barriers.setConstant(T(0.0));
        for (size_t i = 0; i < ev_candidates.size(); ++i) {
            const auto& ev_candidate = ev_candidates[i];
            // a and b are the endpoints of the edge; c is the vertex
            long edge_id = ev_candidate.edge_index;
            int a_id = edges.coeff(edge_id, 0);
            int b_id = edges.coeff(edge_id, 1);
            long c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorX<T> a = vertices_t1.row(a_id);
            Eigen::VectorX<T> b = vertices_t1.row(b_id);
            Eigen::VectorX<T> c = vertices_t1.row(c_id);

            barriers(int(i)) = distance_barrier<T>(a, b, c);
        }
    }
} // namespace opt
} // namespace ccd
