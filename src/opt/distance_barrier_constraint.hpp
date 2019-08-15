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
        typedef DistanceBarrierDiff Diff;

        DistanceBarrierConstraint();
        DistanceBarrierConstraint(const std::string& name);

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        double get_barrier_epsilon () { return m_barrier_epsilon; }
        void set_barrier_epsilon(const double eps) { m_barrier_epsilon = eps; }

        void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXi& group_ids,
            const Eigen::MatrixXd& Uk) override;

        void update_collision_set(const Eigen::MatrixXd& Uk) override;
        void update_active_set(const Eigen::MatrixXd& Uk);

        void compute_constraints(
            const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers);

        void compute_constraints_jacobian(
            const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian);

        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& barriers_hessian);

        template <typename T>
        T distance_barrier(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& c);

        Eigen::VectorXd distance_barrier_grad(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            const Eigen::VectorXd& c);

        Eigen::MatrixXd distance_barrier_hess(const Eigen::VectorXd& a,
            const Eigen::VectorXd& b,
            const Eigen::VectorXd& c);

        const int& number_of_constraints() { return m_num_constraints; }

        const EdgeVertexCandidates& ev_distance_active()
        {
            return m_ev_distance_active;
        }
        size_t ev_distance_active_map(const size_t i)
        {
            return m_ev_distance_active_map[i];
        }
        const std::vector<size_t>& ev_impact_active_map()
        {
            return m_ev_impact_active_map;
        }
        const std::vector<size_t>& ev_inactive_map()
        {
            return m_ev_inactive_map;
        }
        // Settings
        // ----------
        /// @brief initial epsilon to use in barrier function
        double custom_inital_epsilon;
        /// @brief active constraints have distances < scale * barrier_epsilon
        double active_constraint_scale;

    protected:
        /// @brief current number of constraints (includes innactive ones)
        int m_num_constraints;
        /// @brief current barrier epsilon
        double m_barrier_epsilon;

        /// @brief edge-vertex candidates for Impact and Distance evaluations
        EdgeVertexCandidates m_ev_candidates;

        /// @brief edge-vertex candidates with an active barrier
        ///  distance < m_barrier_epsilon * active_cstr_scale
        EdgeVertexCandidates m_ev_distance_active;

        /// @brief map from m_ev_distance_active to m_ev_candidate
        std::vector<size_t> m_ev_distance_active_map;

        /// @brief map from ev_impacts to m_ev_candidate
        std::vector<size_t> m_ev_impact_active_map;

        /// @brief list of inactive constraints
        std::vector<size_t> m_ev_inactive_map;
    };

    template <typename T>
    T point_to_edge_sq_distance(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c);

} // namespace opt
} // namespace ccd

#include "distance_barrier_constraint.tpp"
