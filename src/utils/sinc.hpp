#pragma once

#include <Eigen/Geometry>

#include <autodiff/autodiff.h>
#include <interval/interval.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

double sinc(const double& x);

Interval sinc(const Interval& x);

/// Compute the L2 norm (∑|xᵢ|²)
template <typename T> T norm(const VectorMax3<T>& x)
{
    // Do an explicit abs to avoid possible problems with intervals
    VectorMax3<T> absx(x.size());
    for (int i = 0; i < x.size(); i++) {
        absx(i) = abs(x(i));
    }
    return sqrt(absx.dot(absx));
}

/// Compute sinc(||x||) for x ∈ Rⁿ
template <typename T> T sinc_normx(const VectorMax3<T>& x)
{
    static_assert(
        !std::is_base_of<DiffScalarBase, T>::value,
        "This version does not work with autodiff!");
    return sinc(norm(x));
}

/// Compute ∇sinc(||x||) for x ∈ Rⁿ
VectorMax3d sinc_normx_grad(const VectorMax3d& x);

/// Compute ∇²sinc(||x||) for x ∈ Rⁿ
MatrixMax3d sinc_normx_hess(const VectorMax3d& x);

template <typename Scalar, typename Gradient>
DScalar1<Scalar, Gradient>
sinc_normx(const VectorMax3<DScalar1<Scalar, Gradient>>& x)
{
    const int m = DiffScalarBase::getVariableCount(), n = x.size();

    // Extract the vector of values
    VectorMax3<Scalar> x_vals(n);
    for (int i = 0; i < n; i++) {
        x_vals(i) = x(i).getValue();
    }

    // Compute the value and gradient (without chain-rule)
    Scalar value = sinc_normx(x_vals);
    VectorMax3d grad = sinc_normx_grad(x_vals);

    // Apply chain-rule
    Gradient full_grad = Gradient::Zero(m);
    for (int j = 0; j < n; j++) {
        full_grad += grad(j) * x(j).getGradient();
    }

    return DScalar1<Scalar, Gradient>(value, full_grad);
}

template <typename Scalar, typename Gradient, typename Hessian>
DScalar2<Scalar, Gradient, Hessian>
sinc_normx(const VectorMax3<DScalar2<Scalar, Gradient, Hessian>>& x)
{
    const int m = DiffScalarBase::getVariableCount(), n = x.size();

    // Extract the vector of values
    VectorMax3<Scalar> x_vals(n);
    for (int i = 0; i < n; i++) {
        x_vals(i) = x(i).getValue();
    }

    // Compute the value, gradient, and hessian (without chain-rule)
    Scalar value = sinc_normx(x_vals);
    VectorMax3d grad = sinc_normx_grad(x_vals);
    MatrixMax3d hess = sinc_normx_hess(x_vals);

    // Apply chain-rule
    Gradient full_grad = Gradient::Zero(m);
    for (int k = 0; k < n; k++) {
        full_grad += grad(k) * x(k).getGradient();
    }

    Hessian full_hess = Hessian::Zero(m, m);
    for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
            full_hess += hess(k, l) * x(l).getGradient()
                * x(k).getGradient().transpose();
        }
        full_hess += grad(k) * x(k).getHessian();
    }

    return DScalar2<Scalar, Gradient, Hessian>(value, full_grad, full_hess);
}

} // namespace ipc::rigid
