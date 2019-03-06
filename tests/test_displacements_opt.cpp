#include <iostream>

#include <Eigen/Core>
#include <catch.hpp>
#include <nlopt.hpp>

#include <opt/displacements_nlopt.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Displacment optimization test", "[opt]")
{
    Eigen::MatrixX2d V;
    V.resize(4, 2);
    Eigen::MatrixX2d U;
    U.resize(4, 2);
    Eigen::MatrixX2i E;
    E.resize(2, 2);
    double volume_epsilon = 1.0;
    DetectionMethod ccd_detection_method = DetectionMethod::BRUTE_FORCE;
    Eigen::MatrixX2d Uopt;
    Uopt.resize(4, 2);
    Eigen::MatrixX2d actual_Uopt;
    actual_Uopt.resize(4, 2);

    SECTION("Horizontal edge falling")
    {
        V << -1, 1, 1, 1, -2, 0, 2, 0;
        U << 0, -2, 0, -2, 0, 0, 0, 0;
        E << 0, 1, 2, 3;
        displacements_nlopt_step(
            V, U, E, volume_epsilon, ccd_detection_method, Uopt);
        actual_Uopt << 0, -1, 0, -1, 0, 0, 0, 0;
        CHECK((Uopt - actual_Uopt).squaredNorm() == Approx(0.0).margin(1e-8));
    }
    // SECTION("Vertical edge falling")
    // {
    //     V << 0, 2, 0, 1, -1, 0, 1, 0;
    //     U << 0, -3, 0, -3, 0, 0, 0, 0;
    //     E << 0, 1, 2, 3;
    //     displacements_nlopt_step(
    //         V, U, E, volume_epsilon, ccd_detection_method, Uopt);
    //     V(2, 0) *= 2;
    //     V(3, 0) *= 2;
    //     displacements_nlopt_step(
    //         V, U, E, volume_epsilon, ccd_detection_method, actual_Uopt);
    //     CHECK((Uopt - actual_Uopt).squaredNorm() ==
    //     Approx(0.0).margin(1e-8));
    // }
}
