#pragma once

#include <Eigen/Core>

#include <opt/collision_constraint.hpp>

#include <ccd/collision_detection.hpp>

namespace ccd {
namespace opt {

    class BarrierConstraint : public CollisionConstraint {
    public:
        BarrierConstraint();

        void resetBarrierEpsilon();
        void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk) override;

        void compute_constraints(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& barriers) override;

        void compute_constraints_jacobian(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            Eigen::MatrixXd& barriers_jacobian) override;

        void compute_constraints_hessian(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            std::vector<Eigen::MatrixXd>& barriers_hessian) override;

        // Settings
        // ----------
        double barrier_epsilon;
    };

    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::VectorXd& volumes);

    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& volumes);

    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& volumes_hess);

} // namespace opt
} // namespace ccd
