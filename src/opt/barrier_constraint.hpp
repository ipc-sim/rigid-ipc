#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {
namespace opt {

    class BarrierConstraint {
    public:
        BarrierConstraint();

        void compute_constraints(const Eigen::MatrixX2d vertices,
            const Eigen::MatrixX2i edges, const Eigen::MatrixXd& Uk,
            const EdgeEdgeImpacts& ee_impacts, Eigen::VectorXd& volumes);

        void compute_constraints_jacobian(const Eigen::MatrixX2d vertices,
            const Eigen::MatrixX2i edges, const Eigen::MatrixXd& Uk,
            const EdgeEdgeImpacts& ee_impacts, Eigen::MatrixXd& volumes_jac);

        void compute_constraints_hessian(const Eigen::MatrixX2d vertices,
            const Eigen::MatrixX2i edges, const Eigen::MatrixXd& Uk,
            const EdgeEdgeImpacts& ee_impacts,
            std::vector<Eigen::MatrixXd>& volumes_hessian);

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
