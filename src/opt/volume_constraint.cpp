#include "volume_constraint.hpp"

#include <ccd/collision_volume_diff.hpp>
namespace ccd {

namespace opt {

    VolumeConstraint::VolumeConstraint() {}

    void VolumeConstraint::compute_constraints(const Eigen::MatrixX2d vertices,
        const Eigen::MatrixX2i edges, const Eigen::MatrixXd& Uk,
        const EdgeEdgeImpacts& ee_impacts, Eigen::VectorXd& volumes)
    {
        ccd::autodiff::compute_volumes_refresh_toi(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), volume_epsilon, volumes);
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixX2d vertices, const Eigen::MatrixX2i edges,
        const Eigen::MatrixXd& Uk, const EdgeEdgeImpacts& ee_impacts,
        Eigen::MatrixXd& volumes_jac)
    {
        ccd::autodiff::compute_volumes_gradient(vertices, Uk, edges, ee_impacts,
            Eigen::VectorXi(), volume_epsilon, volumes_jac);
    }

    void VolumeConstraint::compute_constraints_hessian(
        const Eigen::MatrixX2d vertices, const Eigen::MatrixX2i edges,
        const Eigen::MatrixXd& Uk, const EdgeEdgeImpacts& ee_impacts,
        std::vector<Eigen::MatrixXd>& volumes_hessian)
    {
        ccd::autodiff::compute_volumes_hessian(vertices, Uk, edges, ee_impacts,
            Eigen::VectorXi(), volume_epsilon, volumes_hessian);
    }

    void compute_volumes_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::VectorXd& volumes)
    {

        auto solver = std::make_unique<VolumeConstraint>();
        solver->volume_epsilon = epsilon;
        return solver->compute_constraints(V, E, U, ee_impacts, volumes);
    }

    void compute_volumes_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& volumes_jac)
    {
        auto solver = std::make_unique<VolumeConstraint>();
        solver->volume_epsilon = epsilon;
        return solver->compute_constraints_jacobian(V, E, U, ee_impacts, volumes_jac);
    }

    void compute_volumes_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& volumes_hess)
    {
        auto solver = std::make_unique<VolumeConstraint>();
        solver->volume_epsilon = epsilon;
        return solver->compute_constraints_hessian(V, E, U, ee_impacts, volumes_hess);
    }

} // namespace opt
} // namespace ccd
