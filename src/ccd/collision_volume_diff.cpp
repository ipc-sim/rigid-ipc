#include "collision_volume_diff.hpp"

#include <autodiff/finitediff.hpp>
#include <ccd/time_of_impact.hpp>
#include <autogen/collision_volume.hpp>

namespace ccd {
namespace autodiff {
    template <typename T>
    T collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        const double epsilon)
    {
        T toi, alpha;
        bool success;
        success = compute_edge_vertex_time_of_impact(Vi, Vj, Vk, Ui, Uj, Uk, toi);
        success = success || temporal_parameterization_to_spatial(Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);

        assert(success);

        return ccd::autogen::collision_volume(Vi, Vj, Ui, Uj, toi, alpha, epsilon);
    }


    Vector8d collision_volume_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const double epsilon)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DScalar(6, 0.0);
        DScalar(7, 0.0);
        DScalar volume(0.0);

        volume = collision_volume(Vi, Vj, Vk, DUi, DUj, DUk, epsilon);

        return volume.getGradient();
    }
    Vector8d collision_volume_grad_fd(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const double epsilon)
    {
        auto f = [&](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);

            double volume = collision_volume(Vi, Vj, Vk, ui, uj, uk, epsilon);
            return volume;
        };

        Vector8d grad;
        Eigen::VectorXd finite_diff;
        Vector8d x;
        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) << 0.0, 0.0;
        ccd::finite_gradient(x, f, finite_diff);
        grad << finite_diff;
        return grad;
    }

    template double collision_volume<double>(
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        const double);

    template DScalar collision_volume<DScalar>(
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        DVector2 const&, DVector2 const&, DVector2 const&,
        const double);

}
}
