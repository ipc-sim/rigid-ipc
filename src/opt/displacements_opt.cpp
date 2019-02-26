#include "displacements_opt.hpp"

#include "barrier.hpp"
#include "newtons_method.hpp"

#include <ccd/collision_detection.hpp>
#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/prune_impacts.hpp>

namespace ccd {
namespace opt {
    void displacements_optimization_step(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const double barrier_s,
        const double barrier_beta, const DetectionMethod ccd_detection_method,
        Eigen::MatrixX2d& Uopt)
    {

        auto unflatten = [](const Eigen::MatrixXd& _x, Eigen::MatrixXd& x) {
            x = _x;
            x.resize(_x.rows() / 2, 2);
        };

        auto barrier = [barrier_s](const double v) -> double {
            return log_barrier(v, barrier_s);
        };

        /// computes the functional f(u) = ||U - Uo||^2 + beta \sum \phi(v)
        /// @param[in] x        : flattened vector with the updated velocities
        /// @returns            : scalar value
        auto f = [&](const Eigen::VectorXd& x) -> double {
            Eigen::MatrixXd u;
            unflatten(x, u);

            EdgeEdgeImpacts ee_impacts;
            Eigen::VectorXi edge_impact_map(E.rows());
            detect_collisions(
                V, u, E, ccd_detection_method, ee_impacts, edge_impact_map);

            Eigen::VectorXd volumes;
            ccd::compute_volumes(
                V, u, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);

            double energy = (U - u).norm();
            double constr = barrier_beta * volumes.unaryExpr(barrier).sum();

            return energy + constr;
        };

        /// computes the gradient of the functional
        /// @param[in] x        : flattened vector with the updated velocities
        /// @returns            : vector representing the gradient
        auto gradient = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::MatrixXd u;
            unflatten(x, u);

            EdgeEdgeImpacts ee_impacts;
            Eigen::VectorXi edge_impact_map(E.rows());
            detect_collisions(
                V, u, E, ccd_detection_method, ee_impacts, edge_impact_map);

            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(V, u, E, ee_impacts,
                edge_impact_map, volume_epsilon, volume_gradient);

            // TODO: apply chain rule instead of this
            Eigen::VectorXd constr_grad;
            constr_grad = volume_gradient.rowwise().sum();

            // TODO: add gradient of energy
            Eigen::VectorXd energy_grad(V.size());
            energy_grad.setZero();
            return constr_grad + energy_grad;
        };

        /// Computes the hessian of the functional
        /// @param[in] x        : flattened vector with the updated velocities
        /// @returns            : Matrix representing the hessian
        auto hessian = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            Eigen::MatrixXd u;
            unflatten(x, u);

            EdgeEdgeImpacts ee_impacts;
            Eigen::VectorXi edge_impact_map(E.rows());
            detect_collisions(
                V, u, E, ccd_detection_method, ee_impacts, edge_impact_map);

            std::vector<Eigen::MatrixXd> volume_hessian;
            ccd::autodiff::compute_volumes_hessian(V, u, E, ee_impacts,
                edge_impact_map, volume_epsilon, volume_hessian);

            // TODO: apply chain rule here
            Eigen::MatrixXd constr_hessian(V.size(), V.size());
            constr_hessian.setZero();

            // TODO: add gradient of energy
            Eigen::MatrixXd energy_hessian(V.size(), V.size());
            energy_hessian.setZero();

            return constr_hessian + energy_hessian;
        };

        /// Evaluates the constraints
        /// @param[in] x        : flattened vector with the updated velocities
        /// @returns            : Matrix representing the hessian
        auto constraint = [&](const Eigen::VectorXd& x) -> bool {
            Eigen::MatrixXd u;
            unflatten(x, u);

            EdgeVertexImpacts ev_impacts;
            EdgeEdgeImpacts ee_impacts;
            ccd::detect_edge_vertex_collisions(
                V, U, E, ev_impacts, ccd_detection_method);
            convert_edge_vertex_to_edge_edge_impacts(E, ev_impacts, ee_impacts);

            return ee_impacts.size() > 0;
        };

        Eigen::VectorXd _x;
        Eigen::MatrixXd u;
        // flatten
        u = U;
        u.resize(U.size(), 1);
        _x = u;

        newtons_method_step(_x, f, gradient, hessian, constraint);
        unflatten(_x, u);
        Uopt = u;
    }

    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map)
    {

        EdgeVertexImpacts ev_impacts;

        ccd::detect_edge_vertex_collisions(
            V, U, E, ev_impacts, detection_method);
        convert_edge_vertex_to_edge_edge_impacts(E, ev_impacts, ee_impacts);
        ccd::prune_impacts(ee_impacts, edge_impact_map);
    }

}

}
