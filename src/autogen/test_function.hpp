#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H

#include <Eigen/Core>

#include <autodiff/autodiff_types.hpp>

namespace ccd {
namespace autogen {

    template<typename T>
    T _test_function_(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Matrix<T, 2, 1>& Ui,
        const Eigen::Matrix<T, 2, 1>& Uj,
        const Eigen::Matrix<T, 2, 1>& Uk,
        const Eigen::Matrix<T, 2, 1>& Ul,
        const double epsilon);

    double test_function(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon);

    Vector8d test_function_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon);
}
}
#endif
