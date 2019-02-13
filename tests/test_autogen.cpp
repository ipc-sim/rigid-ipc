#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <autogen/test_function.hpp>

using namespace ccd;

TEST_CASE("Test Autogen Test", "[autogen][dummy]")
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

    SECTION("One displacement different than zero -> nonzero energy")
    {
        epsilon = 1.0;
        Ui.setConstant(1.0);
        double e = ccd::autogen::test_function(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

        REQUIRE(e == 2.0);
    }
    SECTION("One position different than zero -> zero energy")
    {
        epsilon = 1.0;
        Xi.setConstant(1.0);
        double e = ccd::autogen::test_function(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

        REQUIRE(e == 0.0);
    }
    SECTION("Energy square of input")
    {
        epsilon = 1.0;
        Ui.setConstant(2.0);
        double e = ccd::autogen::test_function(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

        REQUIRE(e == 8.0);
    }
    SECTION("One displacement different than zero -> one entry in gradient energy")
    {
        epsilon = 1.0;
        Vector8d expected;
        Ui(0) = 1.0;
        expected << 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        Vector8d actual = ccd::autogen::test_function_grad(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

        REQUIRE(actual.rows() == 8);
        REQUIRE(!std::isnan(actual.squaredNorm()));
        REQUIRE(compare_gradient(actual, expected));
    }
    SECTION("One displacement different than zero -> one entry in gradient energy")
    {
        epsilon = 1.0;
        Vector8d expected;
        for (long j = 0; j < 8; j++) {
            for (size_t i = 0; i < 4; i++) {
                X[i]->setZero();
                U[i]->setZero();
            }
            expected.setZero();
            (*U[int(j / 2)])[j % 2] = 1.0;

            expected.row(j).setConstant(2.0);
            Vector8d actual = ccd::autogen::test_function_grad(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);

            REQUIRE(actual.rows() == 8);
            REQUIRE(!std::isnan(actual.squaredNorm()));
            REQUIRE(compare_gradient(actual, expected));
        }
    }
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
            return ccd::autogen::test_function(Xi, Xj, Xk, Xl, ui, uj, uk, ul, epsilon);
        };

        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) = Ul;

        // compare finite differences
        finite_gradient(x, f, expected);
        Vector8d actual = ccd::autogen::test_function_grad(Xi, Xj, Xk, Xl, Ui, Uj, Uk, Ul, epsilon);
        REQUIRE(expected.rows() == 8);
        REQUIRE(actual.rows() == 8);
        REQUIRE(!std::isnan(actual.squaredNorm()));
        REQUIRE(compare_gradient(actual, expected));

    }
}
