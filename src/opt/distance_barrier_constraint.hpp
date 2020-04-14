#pragma once

#include <Eigen/Core>

#include <opt/collision_constraint.hpp>

#include <autodiff/autodiff_types.hpp>
#include <barrier/barrier.hpp>
#include <ccd/collision_detection.hpp>
#include <ccd/hash_grid.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace opt {

    class DistanceBarrierConstraint : public CollisionConstraint {
    public:
        DistanceBarrierConstraint();
        DistanceBarrierConstraint(const std::string& name);

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;

        virtual void initialize() override;

        double get_barrier_epsilon() const { return m_barrier_epsilon; }
        void set_barrier_epsilon(const double eps) { m_barrier_epsilon = eps; }

        bool has_active_collisions(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<double>& poses_t0,
            const physics::Poses<double>& poses_t1) const;

        void compute_constraints(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<double>& poses,
            Eigen::VectorXd& barriers);

        void construct_active_barrier_set(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<double>& poses,
            Candidates& barriers) const;

        template <typename T>
        void compute_candidates_constraints(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<T>& poses,
            const Candidates& candidates,
            Eigen::Matrix<T, Eigen::Dynamic, 1>& barriers);

        template <typename T>
        T distance_barrier(const T distance, const double eps);

        template <typename T> T distance_barrier(const T distance)
        {
            return distance_barrier(distance, m_barrier_epsilon);
        }

        void compute_distances(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<double>& poses,
            Eigen::VectorXd& distances) const;

        // Settings
        // ----------
        /// @brief initial epsilon to use in barrier function
        double custom_inital_epsilon;

        /// @brief displace barrier evaluation by this value
        double min_distance;

        /// @brief active constraints have distances < scale * barrier_epsilon
        double active_constraint_scale;

    protected:
        BarrierType barrier_type;
        double m_barrier_epsilon;
    };

} // namespace opt
} // namespace ccd

#include "distance_barrier_constraint.tpp"
