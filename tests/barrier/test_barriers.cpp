#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <finitediff.hpp>

#include <barrier/barrier.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Test barriers and their derivatives", "[opt][barrier]")
{
    double s = GENERATE(range(-5, 2));
    s = pow(10, s);

    BarrierType barrier_type =
        GENERATE(BarrierType::IPC, BarrierType::POLY_LOG, BarrierType::SPLINE);

    double x = GENERATE_COPY(take(10, random(s / 2, 0.9 * s))); // ∈ [0, s]
    Vector1d x_vec;
    x_vec << x;

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        x_vec,
        [&](const Eigen::VectorXd& x) {
            return barrier(x[0], s, barrier_type);
        },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier_gradient(x, s, barrier_type);

    CAPTURE(barrier_type, s, x, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));

    fd::finite_gradient(
        x_vec,
        [&](const Eigen::VectorXd& x) {
            return barrier_gradient(x[0], s, barrier_type);
        },
        fgrad);

    grad << barrier_hessian(x, s, barrier_type);

    CAPTURE(barrier_type, s, x, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));
}
