#ifndef AUTODIFF_TYPES_H
#define AUTODIFF_TYPES_H

#include <autodiff/autodiff.h>

namespace ccd {
static constexpr size_t N = 2 * 4;
typedef Eigen::Matrix<double, 8, 1> Vector8d;
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef DScalar2<double, Vector8d, Matrix8d> DScalar;
typedef DScalar::DVector2 DVector2;

static inline DVector2 dvector(size_t index, const Eigen::Matrix<double, 2, 1> &v) {
    return DVector2(DScalar(index, v.x()), DScalar(index + 1, v.y()));
}

}
#endif
