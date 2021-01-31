#include <catch2/catch.hpp>

#include <iostream>
#include <opt/volume_constraint.hpp>

#include <finitediff.hpp>
#include <utils/flatten.hpp>

//-----------------
// Tests
// ---------------------------------------------------

TEST_CASE("Volume Constraint", "[opt][ccd][Volume][!mayfail]")
{
    ccd::opt::VolumeConstraint volume;
    volume.volume_epsilon = 1e-3;
    volume.detection_method = ccd::BRUTE_FORCE;

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Eigen::MatrixX2d displacements(4, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    Eigen::VectorXd v_actual, v_expected;

    SECTION("Vertical Displ Small")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
        v_expected.resize(2);
        v_expected << -0.000166667, -8.33333e-05;
    }

    SECTION("Vertical Displ Long")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        v_expected.resize(4);
        v_expected << -0.0005, -0.00025, 0.0, 0.0;
    }

    SECTION("Horizontal Displ Small")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
        v_expected.resize(4);
        v_expected << -0.000166667, -0.0600001, -0.000166667, -0.0600001;
    }

    SECTION("Horizontal Displ Long")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;
        v_expected.resize(4);
        v_expected << -0.0005, -0.3, -0.0005, -0.3;
    }

    volume.initialize(vertices, edges, Eigen::VectorXi(), displacements);
    volume.compute_constraints(displacements, v_actual);
    CHECK((v_actual - v_expected).squaredNorm() < 1e-10);

    Eigen::MatrixXd jac_actual;
    volume.compute_constraints(displacements, v_actual, jac_actual);
    CHECK((v_actual - v_expected).squaredNorm() < 1e-6);
}

TEST_CASE(
    "Volume Constraint Gradient", "[opt][ccd][Volume][Gradient][!mayfail]")
{
    ccd::opt::VolumeConstraint volume;
    volume.volume_epsilon = 1e-3;
    volume.detection_method = ccd::BRUTE_FORCE;

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Eigen::MatrixX2d displacements(4, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    SECTION("Vertical Displ Small")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
    }

    SECTION("Vertical Displ Long")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;

        displacements.row(2) << 0.0, -1.1;
        displacements.row(3) << 0.0, -1.1;
    }

    SECTION("Horizontal Displ Small")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
    }

    SECTION("Horizontal Displ Long")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;

        displacements.row(2) << 0.0, -1.1;
        displacements.row(3) << 0.0, -1.1;
    }

    volume.initialize(vertices, edges, Eigen::VectorXi(), displacements);

    using namespace ccd::opt;
    Eigen::MatrixXd jac_actual;
    Eigen::VectorXd v_actual;
    volume.compute_constraints(displacements, v_actual, jac_actual);

    Eigen::MatrixXd approx_jac;
    auto f = [&](const Eigen::VectorXd& u) -> Eigen::VectorXd {
        Eigen::VectorXd fx;
        Eigen::MatrixXd x = u;
        ccd::unflatten(x, 2);
        volume.compute_constraints(x, fx);
        return fx;
    };
    Eigen::MatrixXd x = displacements;

    ccd::flatten(x);
    fd::finite_jacobian(x, f, approx_jac);

    REQUIRE(approx_jac.rows() == jac_actual.rows());
    CHECK((approx_jac - jac_actual).norm() / approx_jac.norm() < 1e-6);
}
