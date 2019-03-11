#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

#include <ccd/collision_detection.hpp>
#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

void displacements_optimization_step(const Eigen::MatrixX2d& V,
    const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
    const double volume_epsilon, const double barrier_s,
    const double barrier_beta, const DetectionMethod ccd_detection_method,
    Eigen::MatrixX2d& Uopt)
{

    auto flatten = [](const Eigen::MatrixXd& _x, Eigen::MatrixXd& x) {
        x = _x;
        x.resize(_x.rows() * 2, 1);
    };

    auto unflatten = [](const Eigen::MatrixXd& _x, Eigen::MatrixXd& x) {
        x = _x;
        x.resize(_x.rows() / 2, 2);
    };

    auto barrier = [barrier_s](const double v) -> double {
        return spline_barrier(v, barrier_s);
    };
    auto barrier_gradient = [barrier_s](const double v) -> double {
        return spline_barrier_gradient(v, barrier_s);
    };
    auto barrier_hessian = [barrier_s](const double v) -> double {
        return spline_barrier_hessian(v, barrier_s);
    };

    /// computes the functional f(u) = ||U - Uo||^2 + beta \sum \phi(v)
    /// @param[in] x flattened vector with the updated velocities
    /// @returns scalar value
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

        Eigen::MatrixXd U_flat;
        flatten(U, U_flat);
        double energy = (U_flat - x).squaredNorm();
        double constr = barrier_beta * volumes.unaryExpr(barrier).sum();

        return energy + constr;
    };

    /// computes the gradient of the functional
    /// @param[in] x flattened vector with the updated velocities
    /// @returns vector representing the gradient
    auto gradient = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::MatrixXd u;
        unflatten(x, u);

        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, u, E, ccd_detection_method, ee_impacts, edge_impact_map);

        Eigen::VectorXd volumes;
        ccd::compute_volumes(
            V, u, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);

        Eigen::MatrixXd volume_gradient;
        ccd::autodiff::compute_volumes_gradient(V, u, E, ee_impacts,
            edge_impact_map, volume_epsilon, volume_gradient);

        Eigen::VectorXd constr_grad = barrier_beta
            * (volumes.unaryExpr(barrier_gradient).asDiagonal()
                * volume_gradient.transpose())
                  .colwise()
                  .sum()
                  .transpose();

        Eigen::MatrixXd U_flat;
        flatten(U, U_flat);
        Eigen::VectorXd energy_grad = 2 * (U_flat - x);
        return constr_grad + energy_grad;
    };

    /// Computes the hessian of the functional
    /// @param[in] x flattened vector with the updated velocities
    /// @returns Matrix representing the hessian
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

        Eigen::VectorXd volumes;
        ccd::compute_volumes(
            V, u, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);

        Eigen::MatrixXd volume_gradient;
        ccd::autodiff::compute_volumes_gradient(V, u, E, ee_impacts,
            edge_impact_map, volume_epsilon, volume_gradient);

        Eigen::MatrixXd constr_hessian(V.size(), V.size());
        constr_hessian.setZero();
        for (int edge_index = 0; edge_index < E.rows(); edge_index++) {
            constr_hessian += barrier_hessian(volumes[edge_index])
                    * volume_gradient.col(edge_index)
                    * volume_gradient.col(edge_index).transpose()
                + barrier_gradient(volumes[edge_index])
                    * volume_hessian[size_t(edge_index)];
        }
        constr_hessian *= barrier_beta;

        constr_hessian.diagonal().array() += 2;

        return constr_hessian;
    };

    /// Evaluates the constraints
    /// @param[in] x flattened vector with the updated velocities
    /// @returns Matrix representing the hessian
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

    Eigen::VectorXd _x(U.size());
    _x.setZero();

    newtons_method_step(_x, f, gradient, hessian, constraint);
    Eigen::MatrixXd u;
    unflatten(_x, u);
    Uopt = u;
}
