#include "barrier_constraint.hpp"

#include <ccd/collision_penalty_diff.hpp>

namespace ccd {
namespace opt {

    BarrierConstraint::BarrierConstraint() {}
    BarrierConstraint::BarrierConstraint(const double epsilon)
        : barrier_epsilon(epsilon)
    {
    }

    void BarrierConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk)

    {
        ev_impacts.clear();
        ee_impacts.clear();
        edge_impact_map.resize(edges.rows());
        edge_impact_map.setZero();
        detecteCollisions(vertices, edges, Uk);
        resetBarrierEpsilon();
    }

    void BarrierConstraint::resetBarrierEpsilon()
    {
        // ToDo: Find a better way to provide a starting epsilon
        bool was_barrier_epsilon_reset = false;
        for (long i = 0; i < this->edge_impact_map.size(); i++) {
            if (this->edge_impact_map(i) >= 0) {
                EdgeEdgeImpact impact
                    = this->ee_impacts[this->edge_impact_map(i)];
                double val = impact.time;

                // Penalize the spatial distance too
                // double sum = 0;
                // for (int end = 0; end < 2; end++) {
                //     sum += displacements
                //                .row(edges(impact.impacted_edge_index, end))
                //                .squaredNorm();
                //     sum += displacements
                //                .row(edges(impact.impacting_edge_index, end))
                //                .squaredNorm();
                // }
                // val *= sum;

                if (!was_barrier_epsilon_reset || val < barrier_epsilon) {
                    barrier_epsilon = val;
                    was_barrier_epsilon_reset = true;
                }
            }
        }
        if (!was_barrier_epsilon_reset) {
            barrier_epsilon = 0;
        } else {
            barrier_epsilon += 1e-4;
        }
    }

    void BarrierConstraint::compute_constraints(
        const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        ccd::autodiff::compute_constraints_refresh_toi(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<double>, barriers);
    }

    void BarrierConstraint::compute_constraints_jacobian(
        const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian)
    {
        ccd::autodiff::compute_constraints_gradient(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<DScalar>, barriers_jacobian);
        barriers_jacobian.transposeInPlace();
    }

    void BarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::MatrixXd>& barriers_hessian)
    {
        ccd::autodiff::compute_constraints_hessian(vertices, Uk, edges,
            ee_impacts, Eigen::VectorXi(), barrier_epsilon,
            ccd::autodiff::collision_penalty<DScalar>, barriers_hessian);
    }

    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::VectorXd& barriers)
    {
        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints(V, E, U, barriers);
    }

    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        Eigen::MatrixXd& contraints_jac)
    {

        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints_jacobian(V, E, U, contraints_jac);
    }

    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        std::vector<Eigen::MatrixXd>& contraints_hess)
    {
        auto solver = std::make_unique<BarrierConstraint>();
        solver->barrier_epsilon = epsilon;
        solver->ee_impacts = ee_impacts;
        return solver->compute_constraints_hessian(V, E, U, contraints_hess);
    }

} // namespace opt
} // namespace ccd
