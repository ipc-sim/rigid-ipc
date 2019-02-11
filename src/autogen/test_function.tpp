#include "test_function.hpp"

namespace ccd {
namespace autogen {
    template <typename T>
    T _test_function_(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Matrix<T, 2, 1>& Ui,
        const Eigen::Matrix<T, 2, 1>& Uj,
        const Eigen::Matrix<T, 2, 1>& Uk,
        const Eigen::Matrix<T, 2, 1>& Ul,
        const double epsilon)
    {
        T volume(0.0);

        // {{func_ccode}}

        return volume;
    }

    double test_function(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon)
    {
        return _test_function_(Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, epsilon);
    }

    Vector8d test_function_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DVector2 DUl = dvector(6, Ul);

        DScalar volume(0);

        volume = _test_function_(
            Vi, Vj, Vk, Vl, DUi, DUj, DUk, DUl, epsilon);

        return volume.getGradient();
    }
}
}
