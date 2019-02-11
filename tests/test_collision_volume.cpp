#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <autogen/collision_volume.hpp>

using namespace ccd;

TEST_CASE("Test Autogen Collision Volume", "[autogen][collision_volume]")
{
    Eigen::Vector2d Xi, Xj, Xk, Xl; // positions
    Eigen::Vector2d Ui, Uj, Uk, Ul; // velocities
    double epsilon = 1.0;

    Eigen::Vector2d* X[4] = { &Xi, &Xj, &Xk, &Xl };
    Eigen::Vector2d* U[4] = { &Ui, &Uj, &Uk, &Ul };

    for (size_t i = 0; i < 4; i++) {
        X[i]->setZero();
        U[i]->setZero();
    }

    // NOTE: time of collision and lambda are harcoded

    SECTION("Compare gradient with finite diff")
    {
        srand(0);

        // initialize with random data
        epsilon = 1.0;
        for (size_t i = 0; i < 4; i++) {
            (*X[i]) = Eigen::Vector2d::Random();
            (*U[i]) = Eigen::Vector2d::Random();
        }
        Eigen::VectorXd expected;
        Eigen::VectorXd x(8);

        // wrap our function for easy handling on finite-differences
        auto f = [Xi, Xj, Xk, Xl, epsilon](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);
            Eigen::Vector2d ul = U.segment(6, 2);
            return ccd::autogen::collision_volume(Xi, Xj, Xk, Xl, ui, uj, uk, ul, epsilon);
        };

        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) = Ul;

        // compare finite differences
        finite_gradient(x, f, expected);
        Vector8d actual = ccd::autogen::collision_volume_grad(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

        REQUIRE(expected.rows() == 8);
        REQUIRE(actual.rows() == 8);
        REQUIRE(actual.squaredNorm() != 0.0);
        REQUIRE(!std::isnan(actual.squaredNorm()));
        REQUIRE(compare_gradient(actual, expected));

    }
}
