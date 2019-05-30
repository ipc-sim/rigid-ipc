#include <catch.hpp>

#include <iostream>
#include <opt/volume_constraint.hpp>

//-----------------
// Tests
// ---------------------------------------------------

TEST_CASE("Volume Constraint", "[opt][ccd][Volume]")
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

    Eigen::VectorXd v_actual, v_expected(16);

    SECTION("Vertical Displ Small")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
        v_expected << 0, 0, 0, 0, -0.000166667, -8.33333e-05, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0;
    }

    SECTION("Vertical Displ Long")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        v_expected << 0, 0, 0, 0, -0.0005, -0.00025, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0;
    }

    SECTION("Horizontal Displ Small")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
        v_expected << 0, 0, 0, 0, -0.000166667, -0.0600001, -0.000166667,
            -0.0600001, 0, 0, 0, 0, 0, 0, 0, 0;
    }

    SECTION("Horizontal Displ Long")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;
        v_expected << 0, 0, 0, 0, -0.0005, -0.3, -0.0005, -0.3, 0, 0, 0, 0, 0,
            0, 0, 0;
    }

    volume.initialize(vertices, edges, displacements);
    volume.compute_constraints(displacements, v_actual);
    CHECK((v_actual - v_expected).squaredNorm() < 1e-10);

    Eigen::SparseMatrix<double> jac_actual_sparse;
    Eigen::VectorXi active;
    volume.compute_constraints(
        displacements, v_actual, jac_actual_sparse, active);
     CHECK((v_actual - v_expected).squaredNorm() < 1e-6);

}

TEST_CASE("Volume Constraint Gradient", "[opt][ccd][Volume][Gradient]")
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

    Eigen::MatrixXd jac_actual, jac_expected(16, 8);

    SECTION("Vertical Displ Small")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;
        jac_expected << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000138889, -0.000138889,
            -0, -0, -0.000694444, -0.000694444, 0.00138889, -0, -0, -0, -0, 0,
            -0.000347222, -0.000347222, 0.000833333, -0.000138889, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0;
    }

    SECTION("Vertical Displ Long")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        jac_expected << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00025, -0.00025, -0, -0,
            -0.00025, -0.00025, 0.0005, -0, -0, -0, -0, 0, -0.000125, -0.000125,
            0.0005, -0.00025, -0, 0, 0, 0, -0.0005, -0.0005, 0, 0.001, 0, 0, 0,
            -0, -0.00025, -0.00025, -0, 0.0005, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0;
    }

    SECTION("Horizontal Displ Small")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;

        jac_expected << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000138889, -0.000138889,
            -0, -0, -0.00111111, -0.000277778, 0.00138889, -0, -0, -0,
            0.0833334, -0.0833334, -0.400001, -0.1, 0.600001, 0, 0.000138889,
            -0.000138889, -0, -0, -0.000277778, -0.00111111, -0, 0.00138889, -0,
            -0, 0.0833334, -0.0833334, -0.1, -0.400001, -0, 0.600001, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    }

    SECTION("Horizontal Displ Long")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        jac_expected << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00025, -0.00025, -0, -0,
            -0.0004, -0.0001, 0.0005, -0, -0, -0, 0.25, -0.25, -0.24, -0.06,
            0.6, 0, 0.00025, -0.00025, -0, -0, -0.0001, -0.0004, -0, 0.0005, -0,
            -0, 0.25, -0.25, -0.06, -0.24, -0, 0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0;
    }

    volume.initialize(vertices, edges, displacements);
    volume.compute_constraints_jacobian(displacements, jac_actual);
    CHECK((jac_actual - jac_expected).squaredNorm() < 1e-10);

}

TEST_CASE("Volume Constraint Gradient Sparse", "[opt][ccd][Volume][Gradient]")
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

    Eigen::MatrixXd jac_actual;
    Eigen::SparseMatrix<double> jac_actual_sparse;

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
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

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
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

    }

    volume.initialize(vertices, edges, displacements);
    volume.compute_constraints_jacobian(displacements, jac_actual);
    volume.compute_constraints_jacobian(displacements, jac_actual_sparse);
    CHECK((jac_actual - jac_actual_sparse.toDense()).squaredNorm() < 1e-10);

    Eigen::VectorXd v_actual;
    Eigen::VectorXi active;
    volume.compute_constraints(
        displacements, v_actual, jac_actual_sparse, active);
    CHECK((jac_actual - jac_actual_sparse.toDense()).squaredNorm() < 1e-6);


}

TEST_CASE("Volume Constraint Hessian", "[opt][ccd][Volume][Hessian]")
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

    std::vector<Eigen::SparseMatrix<double>> hess_actual, hess_expected;
    Eigen::MatrixXd hessian(8, 8);
    SECTION("Vertical Displ Small")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;

        hess_expected.resize(2);

        hessian << -1.807e-20, 1.807e-20, 0, 0, -0.000115741, 0.00104167,
            -0.000925926, 0, 1.807e-20, -1.807e-20, 0, 0, -0.00104167,
            0.000115741, 0.000925926, 0, 0, 0, 0, 0, 0.00115741, -0.00115741, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -0.000115741, -0.00104167, 0.00115741, 0,
            -41.6656, -41.6654, -0.00231481, 0, 0.00104167, 0.000115741,
            -0.00115741, 0, -41.6654, -41.6656, -0.00231481, 0, -0.000925926,
            0.000925926, 0, 0, -0.00231481, -0.00231481, 0.00462963, 0, 0, 0, 0,
            0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hessian << 0, 0, 0, 0, -0.000289352, 0.000289352, 0, 0, 0, 0, 0, 0,
            -0.000289352, 0.000289352, 0, 0, 0, 0, -0.000231481, 0.000231481,
            0.000578704, -0.000578704, 0, 0, 0, 0, 0.000231481, -83.3336, 0, 0,
            0, 0, -0.000289352, -0.000289352, 0.000578704, 0, 0.000578704,
            0.000578704, -0.000694444, -0.000462963, 0.000289352, 0.000289352,
            -0.000578704, 0, 0.000578704, 0.000578704, -0.000694444,
            -0.000462963, 0, 0, 0, 0, -0.000694444, -0.000694444, 0.000462963,
            0.000925926, 0, 0, 0, 0, -0.000462963, -0.000462963, 0.000925926,
            -3.61401e-20;
        hess_expected[1] = hessian.sparseView();
    }

    SECTION("Vertical Displ Long")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        hess_expected.resize(4);

        hessian << -2.71051e-20, 2.71051e-20, 0, 0, -0.000125, 0.000125, 0, 0,
            2.71051e-20, -2.71051e-20, 0, 0, -0.000125, 0.000125, 0, 0, 0, 0, 0,
            0, 0.00025, -0.00025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.000125,
            -0.000125, 0.00025, 0, -125, -125, -0.0005, 0, 0.000125, 0.000125,
            -0.00025, 0, -125, -125, -0.0005, 0, 0, 0, 0, 0, -0.0005, -0.0005,
            0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hessian << 0, 0, 0, 0, -6.25e-05, 6.25e-05, 0, 0, 0, 0, 0, 0, -6.25e-05,
            6.25e-05, 0, 0, 0, 0, -0.00025, 0.00025, 0.000125, -0.000125, 0, 0,
            0, 0, 0.00025, -250, 0, 0, 0, 0, -6.25e-05, -6.25e-05, 0.000125, 0,
            0.000125, 0.000125, -0.00025, 0, 6.25e-05, 6.25e-05, -0.000125, 0,
            0.000125, 0.000125, -0.00025, 0, 0, 0, 0, 0, -0.00025, -0.00025,
            0.0005, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[1] = hessian.sparseView();
        hessian << 0, 0, 0, 0, 0, 0.001, 0, -0.001, 0, 0, 0, 0, -0.001, 0, 0,
            0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001, -0.001, 0, 0, 0,
            -0.001, 0, 0.001, 0.0005, 0.0005, 0, -0.001, 0.001, 0, 0, -0.001,
            0.0005, 0.0005, 0, -0.001, 0, 0, 0, 0, 0, 0, 0, 0, -0.001, 0.001, 0,
            0, -0.001, -0.001, 0, 0.002;
        hess_expected[2] = hessian.sparseView();
        hessian << 0, 0, 0, 0, -0.00025, 0.00025, 0, 0, 0, 0, 0, 0, -0.00025,
            0.00025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0005, -0.0005,
            0, 0, -0.00025, -0.00025, 0, 0.0005, 0.00025, 0.00025, 0.0005,
            -0.001, 0.00025, 0.00025, 0, -0.0005, 0.00025, 0.00025, 0.0005,
            -0.001, 0, 0, 0, 0, 0.0005, 0.0005, 0, -0.001, 0, 0, 0, 0, -0.001,
            -0.001, -0.001, 0.003;
        hess_expected[3] = hessian.sparseView();
    }

    SECTION("Horizontal Displ Small")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.6;
        displacements.row(3) << 0.0, -0.6;

        hess_expected.resize(4);

        hessian << -1.807e-20, 1.807e-20, 0, 0, -0.000185185, 0.00111111,
            -0.000925926, 0, 1.807e-20, -1.807e-20, 0, 0, -0.000972222,
            4.62963e-05, 0.000925926, 0, 0, 0, 0, 0, 0.00115741, -0.00115741, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -0.000185185, -0.000972222, 0.00115741,
            0, -106.664, -26.6658, -0.0037037, 0, 0.00111111, 4.62963e-05,
            -0.00115741, 0, -26.6658, -6.6666, -0.000925926, 0, -0.000925926,
            0.000925926, 0, 0, -0.0037037, -0.000925926, 0.00462963, 0, 0, 0, 0,
            0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();

        hessian << 0, 0, 0, 0, -0.333334, 0.333334, 0, 0, 0, 0, 0, 0,
            -0.0833334, 0.0833334, 0, 0, 0, 0, -3.70074e-17, 3.70074e-17,
            0.861112, -0.305556, -0.555556, -0.138889, 0, 0, 3.70074e-17,
            -3.70074e-17, -0.444445, -0.111111, 0.694445, 0, -0.333334,
            -0.0833334, 0.861112, -0.444445, 1.06667, 0.266667, -0.666669, 0,
            0.333334, 0.0833334, -0.305556, -0.111111, 0.266667, 0.0666668,
            -0.166667, 0, 0, 0, -0.555556, 0.694445, -0.666669, -0.166667,
            3.84516e-06, 3.21502e-07, 0, 0, -0.138889, 0, 0, 0, 3.21502e-07,
            -3.21502e-07;
        hess_expected[1] = hessian.sparseView();

        hessian << -1.807e-20, 1.807e-20, 0, 0, -4.62963e-05, 0.000972222, 0,
            -0.000925926, 1.807e-20, -1.807e-20, 0, 0, -0.00111111, 0.000185185,
            0, 0.000925926, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00115741,
            -0.00115741, 0, 0, -4.62963e-05, -0.00111111, 0, 0.00115741,
            -6.6666, -26.6658, 0, -0.000925926, 0.000972222, 0.000185185, 0,
            -0.00115741, -26.6658, -106.664, 0, -0.0037037, 0, 0, 0, 0, 0, 0, 0,
            0, -0.000925926, 0.000925926, 0, 0, -0.000925926, -0.0037037, 0,
            0.00462963;
        hess_expected[2] = hessian.sparseView();

        hessian << 0, 0, 0, 0, -0.0833334, 0.0833334, 0, 0, 0, 0, 0, 0,
            -0.333334, 0.333334, 0, 0, 0, 0, -3.70074e-17, 3.70074e-17,
            0.111111, 0.444445, 0, -0.694445, 0, 0, 3.70074e-17, -3.70074e-17,
            0.305556, -0.861112, 0.138889, 0.555556, -0.0833334, -0.333334,
            0.111111, 0.305556, 0.0666668, 0.266667, 0, -0.166667, 0.0833334,
            0.333334, 0.444445, -0.861112, 0.266667, 1.06667, 0, -0.666669, 0,
            0, 0, 0.138889, 0, 0, -3.21502e-07, 3.21502e-07, 0, 0, -0.694445,
            0.555556, -0.166667, -0.666669, 3.21502e-07, 3.84516e-06;
        hess_expected[3] = hessian.sparseView();
    }

    SECTION("Horizontal Displ Long")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -1.0;
        displacements.row(3) << 0.0, -1.0;

        hess_expected.resize(4);

        hessian << -2.71051e-20, 2.71051e-20, 0, 0, -0.0002, 0.0002, 0, 0,
            2.71051e-20, -2.71051e-20, 0, 0, -5e-05, 5e-05, 0, 0, 0, 0, 0, 0,
            0.00025, -0.00025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0002, -5e-05,
            0.00025, 0, -319.999, -79.9997, -0.0008, 0, 0.0002, 5e-05, -0.00025,
            0, -79.9997, -20.0001, -0.0002, 0, 0, 0, 0, 0, -0.0008, -0.0002,
            0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();

        hessian << 0, 0, 0, 0, -0.12, 0.12, 0, 0, 0, 0, 0, 0, -0.03, 0.03, 0, 0,
            0, 0, 5.55112e-17, -5.55112e-17, 0.15, -0.15, 0, -0.25, 0, 0,
            -5.55112e-17, 5.55112e-17, -2.77556e-17, -6.93889e-18, 0.25, 0,
            -0.12, -0.03, 0.15, -2.77556e-17, 0.384, 0.096, -0.24, 0, 0.12,
            0.03, -0.15, -6.93889e-18, 0.096, 0.024, -0.0600001, 0, 0, 0, 0,
            0.25, -0.24, -0.0600001, 9.16669e-08, 2.08333e-07, 0, 0, -0.25, 0,
            0, 0, 2.08333e-07, -2.08333e-07;
        hess_expected[1] = hessian.sparseView();

        hessian << -2.71051e-20, 2.71051e-20, 0, 0, -5e-05, 5e-05, 0, 0,
            2.71051e-20, -2.71051e-20, 0, 0, -0.0002, 0.0002, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0.00025, -0.00025, 0, 0, -5e-05, -0.0002, 0,
            0.00025, -20.0001, -79.9997, 0, -0.0002, 5e-05, 0.0002, 0, -0.00025,
            -79.9997, -319.999, 0, -0.0008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -0.0002, -0.0008, 0, 0.001;
        hess_expected[2] = hessian.sparseView();

        hessian << 0, 0, 0, 0, -0.03, 0.03, 0, 0, 0, 0, 0, 0, -0.12, 0.12, 0, 0,
            0, 0, 5.55112e-17, -5.55112e-17, 6.93889e-18, 2.77556e-17, 0, -0.25,
            0, 0, -5.55112e-17, 5.55112e-17, 0.15, -0.15, 0.25, 0, -0.03, -0.12,
            6.93889e-18, 0.15, 0.024, 0.096, 0, -0.0600001, 0.03, 0.12,
            2.77556e-17, -0.15, 0.096, 0.384, 0, -0.24, 0, 0, 0, 0.25, 0, 0,
            -2.08333e-07, 2.08333e-07, 0, 0, -0.25, 0, -0.0600001, -0.24,
            2.08333e-07, 9.16669e-08;
        hess_expected[3] = hessian.sparseView();
    }

    volume.initialize(vertices, edges, displacements);
    volume.compute_constraints_hessian(displacements, hess_actual);

    CHECK(hess_actual.size() == hess_expected.size());
    for (size_t i = 0; i < hess_actual.size(); i++) {
        CHECK((hess_actual[i].toDense() - hess_expected[i].toDense()).squaredNorm()
            < 1e-6);
    }

}
