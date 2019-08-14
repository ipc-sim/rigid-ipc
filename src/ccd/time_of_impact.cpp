#include "time_of_impact.hpp"

#include <autodiff/finitediff.hpp>

namespace ccd {
namespace autodiff {


    void compute_edge_vertex_time_of_impact_grad(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk, Vector8d& grad)
    {

        // All definitions using DScalar must be done after setVariableCount
        // NOTE: gradient is computed over 4 vertex (2 for each edge)
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DScalar(6, 0.0);
        DScalar(7, 0.0);

        DScalar toi(0.0);

        compute_edge_vertex_time_of_impact(Vi, Vj, Vk, DUi, DUj, DUk, toi);
        grad = toi.getGradient();
    }

    void compute_edge_vertex_time_of_impact_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk, Vector8d& grad)
    {
        auto f = [&](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);

            double toi;
            bool success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vk, ui, uj, uk, toi);
            return success ? toi : 0.0;
        };

        Eigen::VectorXd finite_diff;
        Vector8d x;
        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) << 0.0, 0.0;
        ccd::finite_gradient(x, f, finite_diff);
        grad << finite_diff;
    }

    template bool compute_edge_vertex_time_of_impact<double>(
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, double&);

    template bool temporal_parameterization_to_spatial<double>(
        const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Vector2T<double>& Ui,
        const Vector2T<double>& Uj, const Vector2T<double>& Uk,
        const double& toi, double& alpha);
} // namespace autodiff

} // namespace ccd
