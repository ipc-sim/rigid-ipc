#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Core>
#include <catch2/catch.hpp>
#include <nlopt.hpp>

/**
 * @brief \f$\min \frac{1}{2}x^Tx\f$
 */
double test_func(const std::vector<double>& x, std::vector<double>& grad,
    void* /* my_func_data */)
{
    Eigen::VectorXd X
        = Eigen::Map<const Eigen::VectorXd>(x.data(), long(x.size()));
    if (!grad.empty()) {
        Eigen::VectorXd gradX = X;
        assert(size_t(gradX.size()) == size_t(grad.size()));
        Eigen::VectorXd::Map(&grad[0], grad.size()) = gradX;
    }
    return X.squaredNorm() / 2;
}

TEST_CASE("Simple tests of NLOPT", "[opt][nlopt]")
{
    const int DIM = 10;
    nlopt::algorithm algo = GENERATE(nlopt::LD_MMA, nlopt::LD_SLSQP);
    nlopt::opt opt(algo, DIM);
    std::vector<double> lb(DIM, -HUGE_VAL); // lower bounds
    lb[0] = GENERATE(0, 1);
    opt.set_lower_bounds(lb);
    opt.set_min_objective(test_func, nullptr);
    opt.set_xtol_rel(1e-8);
    opt.set_xtol_abs(1e-8);
    std::vector<double> x(
        DIM, GENERATE(1, 10, 100, 1000)); // some initial guess
    double minf; // the minimum objective value, upon return
    nlopt::result result = opt.optimize(x, minf);
    CHECK(x[0] == Approx(lb[0]).margin(1e-8));
    CHECK(minf == Approx(lb[0] * lb[0] / 2).margin(1e-8));
}
