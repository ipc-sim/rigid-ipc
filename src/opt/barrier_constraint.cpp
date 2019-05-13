#include "barrier_constraint.hpp"

#include <ccd/collision_penalty_diff.hpp>

namespace ccd {
namespace opt {

    BarrierConstraint::BarrierConstraint() {}

    void BarrierConstraint::compute_constraints(const Eigen::MatrixX2d vertices,
        const Eigen::MatrixX2i edges, const Eigen::MatrixXd& Uk,
        const EdgeEdgeImpacts& ee_impacts, Eigen::VectorXd& constraints)
    {
        ccd::autodiff::compute_constraints_refresh_toi(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<double>, constraints);
    }

    void BarrierConstraint::compute_constraints_jacobian(
        const Eigen::MatrixX2d vertices, const Eigen::MatrixX2i edges,
        const Eigen::MatrixXd& Uk, const EdgeEdgeImpacts& ee_impacts,
        Eigen::MatrixXd& constraints_jac)
    {
        ccd::autodiff::compute_constraints_gradient(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<DScalar>, constraints_jac);
    }
    void BarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixX2d vertices, const Eigen::MatrixX2i edges,
        const Eigen::MatrixXd& Uk, const EdgeEdgeImpacts& ee_impacts,
        std::vector<Eigen::MatrixXd>& constraints_hessian)
    {
        ccd::autodiff::compute_constraints_hessian(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<DScalar>, constraints_hessian);
    }

    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::VectorXd& constraints)
    {
        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        return solver->compute_constraints(V, E, U, ee_impacts, constraints);
    }

    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& contraints_jac)
    {

        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        return solver->compute_constraints_jacobian(
            V, E, U, ee_impacts, contraints_jac);
    }

    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& contraints_hess)
    {
        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        return solver->compute_constraints_hessian(
            V, E, U, ee_impacts, contraints_hess);
    }
} // namespace opt
} // namespace ccd
