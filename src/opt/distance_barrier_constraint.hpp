#pragma once

#include <Eigen/Core>

#include <opt/collision_constraint.hpp>

#include <autodiff/autodiff_types.hpp>
#include <ccd/collision_detection.hpp>
#include <ccd/hash_grid.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    class DistanceBarrierConstraint : public CollisionConstraint {
    public:
        typedef AutodiffType<6> Diff;
        enum CollisionCheck { HAS_COLLISION = 0, NO_COLLISIONS };

        DistanceBarrierConstraint();
        DistanceBarrierConstraint(const std::string& name);

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        double get_barrier_epsilon() { return m_barrier_epsilon; }
        void set_barrier_epsilon(const double eps) { m_barrier_epsilon = eps; }

        EdgeVertexImpacts initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXi& group_ids,
            const Eigen::MatrixXd& Uk) override;

        CollisionCheck get_active_barrier_set(
            const Eigen::MatrixXd& Uk, EdgeVertexCandidates& ev_barriers);

        void compute_constraints(
            const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers);

        void compute_constraints_jacobian(
            const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian);

        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& barriers_hessian);

        template <typename T>
        void compute_candidates_constraints(
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Uk,
            const EdgeVertexCandidates& ev_candidates,
            Eigen::Matrix<T, Eigen::Dynamic, 1>& barriers);

        void compute_candidates_constraints_jacobian(const Eigen::MatrixXd& Uk,
            const EdgeVertexCandidates& ev_candidates,
            Eigen::MatrixXd& barriers_jacobian);

        void compute_candidates_constraints_hessian(const Eigen::MatrixXd& Uk,
            const EdgeVertexCandidates& ev_candidates,
            std::vector<Eigen::SparseMatrix<double>>& barriers_hessian);

        template <typename T>
        T distance_barrier(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& c);

        template <typename T>
        T distance_barrier(const T distance);

        double distance_barrier_grad(const double distance);

        Eigen::VectorXd distance_barrier_grad(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            const Eigen::VectorXd& c);

        Eigen::MatrixXd distance_barrier_hess(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            const Eigen::VectorXd& c);

        // Settings
        // ----------
        /// @brief initial epsilon to use in barrier function
        double custom_inital_epsilon;
        /// @brief active constraints have distances < scale * barrier_epsilon
        double active_constraint_scale;

    protected:
        double m_barrier_epsilon;
    };

    template <typename T>
    T point_to_edge_sq_distance(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c);

} // namespace opt
} // namespace ccd

#include "distance_barrier_constraint.tpp"
