#include "volume_constraint.hpp"

#include <ccd/collision_volume_diff.hpp>
namespace ccd {

namespace opt {

    VolumeConstraint::VolumeConstraint()
        : volume_epsilon(1E-3)
    {
    }

    void VolumeConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk)

    {
        ev_impacts.clear();
        ee_impacts.clear();
        edge_impact_map.resize(edges.rows());
        edge_impact_map.setZero();
        detecteCollisions(vertices, edges, Uk);
    }

    void VolumeConstraint::compute_constraints(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& volumes)
    {
        ccd::autodiff::compute_volumes_refresh_toi(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), volume_epsilon, volumes);
    }

    void VolumeConstraint::compute_constraints_jacobian(
        const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& volumes_jac)
    {
        ccd::autodiff::compute_volumes_gradient(vertices, Uk, edges, ee_impacts,
            Eigen::VectorXi(), volume_epsilon, volumes_jac);
        volumes_jac.transposeInPlace();
    }

    void VolumeConstraint::compute_constraints_hessian(
        const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk,
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
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints(V, E, U, volumes);
    }

    void compute_volumes_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& volumes_jac)
    {
        auto solver = std::make_unique<VolumeConstraint>();
        solver->volume_epsilon = epsilon;
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints_jacobian(V, E, U, volumes_jac);
    }

    void compute_volumes_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& volumes_hess)
    {
        auto solver = std::make_unique<VolumeConstraint>();
        solver->volume_epsilon = epsilon;
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints_hessian(V, E, U, volumes_hess);
    }

} // namespace opt
} // namespace ccd
