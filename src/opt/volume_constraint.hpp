#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/collision_constraint.hpp>

namespace ccd {

namespace opt {

    class VolumeConstraint : public CollisionConstraint {
    public:
        VolumeConstraint();
        void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk) override;

        void compute_constraints(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& volumes) override;

        void compute_constraints_jacobian(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            Eigen::MatrixXd& volumes_jac) override;

        void compute_constraints_hessian(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            std::vector<Eigen::MatrixXd>& volumes_hessian) override;

        // Settings
        // ----------
        double volume_epsilon;
    };

    void compute_volumes_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::VectorXd& volumes);

    void compute_volumes_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& volumes);

    void compute_volumes_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& volumes_hess);

} // namespace opt
} // namespace ccd
