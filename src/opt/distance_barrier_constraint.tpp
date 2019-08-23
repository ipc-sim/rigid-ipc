#pragma once

#include "distance_barrier_constraint.hpp"
#include <opt/barrier.hpp>

namespace ccd {
namespace opt {

    template <typename T>
    T point_to_edge_sq_distance(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> ab = b - a;
        Eigen::Matrix<T, Eigen::Dynamic, 1> ac = c - a;
        Eigen::Matrix<T, Eigen::Dynamic, 1> bc = c - b;

        T e = ac.dot(ab);

        // Handle cases where c projects outside ab
        if (e <= T(0.0)) {
            return ac.dot(ac);
        }
        T f = ab.dot(ab);
        if (e >= f) {
            return bc.dot(bc);
        }

        // Handle cases where c projects onto ab
        T g = ac.dot(ac);
        return g - e * e / f;
    }

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c)
    {
        T distance = sqrt(point_to_edge_sq_distance<T>(a, b, c));
        return opt::spline_barrier<T>(distance, m_barrier_epsilon);
    }

    template <typename T>
    void DistanceBarrierConstraint::compute_candidates_constraints(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Uk,
        const EdgeVertexCandidates& ev_candidates,
        Eigen::Matrix<T, Eigen::Dynamic, 1>& barriers)
    {
        typedef const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
        typedef const Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXd;
        // distance barrier is evaluated at end-positions
        MatrixXd vertices_t1 = vertices.cast<T>() + Uk;

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
            VectorXd a = vertices_t1.row(a_id);
            VectorXd b = vertices_t1.row(b_id);
            VectorXd c = vertices_t1.row(c_id);

            barriers(int(i)) = distance_barrier<T>(a, b, c);
        }
    }
} // namespace opt
} // namespace ccd
