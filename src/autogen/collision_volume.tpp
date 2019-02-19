#include "collision_volume.hpp"

#include <iostream>

#include <autodiff/finitediff.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    void time_of_impact_coeff(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        T& a, T& b, T& c)
    {
        // {{toi_abc_ccode}}
    }

    template <typename T>
    T _collision_volume_(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const double epsilon)
    {
        using namespace std;
        double alpha = 0.5, toi = 0.5;
        T volume(0.0);

        // {{volume_ccode}}

        return volume;
    }

    double collision_volume(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const double epsilon)
    {
        return _collision_volume_(Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, epsilon);
    }

    Vector8d collision_volume_grad(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const double epsilon)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DVector2 DUl = dvector(6, Ul);

        DScalar volume(0.0);

        volume
            = _collision_volume_(Vi, Vj, Vk, Vl, DUi, DUj, DUk, DUl, epsilon);

        return volume.getGradient();
    }

    template void time_of_impact_coeff<double>(
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        double&, double&, double&);

    template void time_of_impact_coeff<DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> > >(
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >, 2, 1, 0, 2, 1> const&,
        DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >&,
        DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >&,
        DScalar1<double, Eigen::Matrix<double, 8, 1, 0, 8, 1> >&);


}
}
