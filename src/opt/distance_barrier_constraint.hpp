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

        bool is_barrier() override { return true; }
        double get_barrier_epsilon() override { return barrier_epsilon; }
        void set_barrier_epsilon(const double eps) override
        {
            barrier_epsilon = eps;
        }

        void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::MatrixXd& Uk) override;

        int number_of_constraints() override;

        virtual void detectCollisions(const Eigen::MatrixXd& Uk) override;

        void compute_constraints(
            const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers) override;

        void compute_constraints_jacobian(const Eigen::MatrixXd& Uk,
            Eigen::MatrixXd& barriers_jacobian) override;

        void compute_constraints_hessian(const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
            override;

        void compute_constraints_and_derivatives(const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk,
            Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

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

        // Settings
        // ----------
        double custom_inital_epsilon;
        bool use_hash_grid;

    protected:
        double barrier_epsilon;
        int num_active_constraints;
        EdgeVertexCandidates ev_distance_candidates;

        Eigen::VectorXi constraint_map; ///< map total number of possible
                                        ///< constraints to active set

        void compute_candidate_intersections_hashgrid(
            const Eigen::MatrixXd& vertices_t1,
            EdgeVertexCandidates& candidates) const;

        void compute_candidate_intersections_brute_force(
            const Eigen::MatrixXd& vertices_t1,
            EdgeVertexCandidates& candidates) const;
    };

    template <typename T>
    T point_to_edge_sq_distance(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c);

} // namespace opt
} // namespace ccd
