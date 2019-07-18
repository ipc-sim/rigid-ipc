#include <catch2/catch.hpp>

#include <iostream>
#include <opt/time_barrier_constraint.hpp>

//-----------------
// Tests
//-----------------

TEST_CASE("Time Barrier Constraint", "[opt][ccd][Barrier]")
{
    ccd::opt::TimeBarrierConstraint barrier;
    barrier.detection_method = ccd::BRUTE_FORCE;

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Eigen::MatrixX2d displacements(4, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    Eigen::VectorXd v_actual, v_expected(2);

    SECTION("Vertical Displ Small Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;

        v_expected << 82.6676, 41.3338;
    }

    SECTION("Vertical Displ Large Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;
        v_expected << 1.04918, 0.52459;
    }

    SECTION("Horizontal Displ Small Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;
        v_expected.resize(4);
        v_expected << 82.6676, 49.6005, 82.6676, 49.6005;
    }

    SECTION("Horizontal Displ Big Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;
        v_expected.resize(4);
        v_expected << 1.04918, 0.629508, 1.04918, 0.629508;
    }

    REQUIRE(
        barrier.initial_epsilon == ccd::opt::InitialBarrierEpsilon::MIN_TOI);
    barrier.initialize(vertices, edges, displacements * 2.0);

    barrier.compute_constraints(displacements, v_actual);
    CHECK((v_actual - v_expected).squaredNorm() < 1e-6);

    Eigen::MatrixXd jac_actual;
    std::vector<Eigen::SparseMatrix<double>> hess_actual;
    barrier.compute_constraints_and_derivatives(
        displacements, v_actual, jac_actual, hess_actual);
    CHECK((v_actual - v_expected).squaredNorm() < 1e-6);

    //    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
    //    Eigen::DontAlignCols,
    //        ", ", ", ", "", "", " << ", ";");
    //    std::cout << v_actual.format(CommaInitFmt) << std::endl;
    //    std::cout << v_expected.format(CommaInitFmt) << std::endl;
}

TEST_CASE("Time Barrier Constraint Gradient", "[opt][ccd][Barrier][Gradient]")
{
    ccd::opt::TimeBarrierConstraint barrier;
    barrier.detection_method = ccd::BRUTE_FORCE;

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Eigen::MatrixX2d displacements(4, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    Eigen::MatrixXd jac_actual, jac_expected(2, 8);

    SECTION("Vertical Displ Small Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;
        jac_expected << -0, -0, -0, -0, 41749.71988, 41749.71988, -83499.43976,
            -0, -0, -0, -0, -0, 20874.85994, 20874.85994, -41749.71988, -0;
    }

    SECTION("Vertical Displ Large Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;

        jac_expected << -0, -0, -0, -0, 17.9163, 17.9163, -35.8327, -0, -0, -0,
            -0, -0, 8.95817, 8.95817, -17.9163, -0;
    }

    SECTION("Horizontal Displ Small Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;

        jac_expected.resize(4, 8);
        jac_expected << -0, -0, -0, -0, 66799.55181, 16699.88795, -83499.43976,
            -0, -0, -0, -0, -0, 40079.73109, 10019.93277, -50099.66386, -0, -0,
            -0, -0, -0, 16699.88795, 66799.55181, -0, -83499.43976, -0, -0, -0,
            -0, 10019.93277, 40079.73109, -0, -50099.66386;
    }

    SECTION("Horizontal Displ Big Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;

        jac_expected.resize(4, 8);
        jac_expected << -0, -0, -0, -0, 28.6661, 7.16653, -35.8327, -0, -0, -0,
            -0, -0, 17.1997, 4.29992, -21.4996, -0, -0, -0, -0, -0, 7.16653,
            28.6661, -0, -35.8327, -0, -0, -0, -0, 4.29992, 17.1997, -0,
            -21.4996;
    }

    REQUIRE(
        barrier.initial_epsilon == ccd::opt::InitialBarrierEpsilon::MIN_TOI);
    barrier.initialize(vertices, edges, displacements * 2.0);

    barrier.compute_constraints_jacobian(displacements, jac_actual);
    CHECK((jac_actual - jac_expected).squaredNorm() < 1e-6);

    Eigen::VectorXd v_actual;
    std::vector<Eigen::SparseMatrix<double>> hess_actual;
    barrier.compute_constraints_and_derivatives(
        displacements, v_actual, jac_actual, hess_actual);
    CHECK((jac_actual - jac_expected).squaredNorm() < 1e-6);

    //    Eigen::IOFormat CommaInitFmt(
    //        10, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
    //    std::cout << v_actual.format(CommaInitFmt) << std::endl;
    //    std::cout << v_expected.format(CommaInitFmt) << std::endl;
}

TEST_CASE("Time Barrier Constraint Hessian", "[opt][ccd][Barrier][Hessian]")
{
    ccd::opt::TimeBarrierConstraint barrier;
    barrier.detection_method = ccd::BRUTE_FORCE;

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

    SECTION("Vertical Displ Small Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;

        hess_expected.resize(2);
        hessian << 0, 0, 0, 0, 41833.3866539627, -41833.3866539627, 0, 0, 0, 0,
            0, 0, 41833.3866539627, -41833.3866539627, 0, 0, 0, 0, 0, 0,
            -83666.7733079253, 83666.7733079253, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            41833.3866539627, 41833.3866539627, -83666.7733079253, 0,
            41750168.7927318, 41750168.7927318, -83500337.5854636, 0,
            -41833.3866539627, -41833.3866539627, 83666.7733079253, 0,
            41750168.7927318, 41750168.7927318, -83500337.5854636, 0, 0, 0, 0,
            0, -83500337.5854636, -83500337.5854636, 167000675.170927, 0, 0, 0,
            0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hess_expected.resize(2);
        hessian << 0, 0, 0, 0, 20916.6933269813, -20916.6933269813, 0, 0, 0, 0,
            0, 0, 20916.6933269813, -20916.6933269813, 0, 0, 0, 0, 0, 0,
            -41833.3866539627, 41833.3866539627, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            20916.6933269813, 20916.6933269813, -41833.3866539627, 0,
            20875084.3963659, 20875084.3963659, -41750168.7927318, 0,
            -20916.6933269813, -20916.6933269813, 41833.3866539627, 0,
            20875084.3963659, 20875084.3963659, -41750168.7927318, 0, 0, 0, 0,
            0, -41750168.7927318, -41750168.7927318, 83500337.5854636, 0, 0, 0,
            0, 0, 0, 0, 0, 0;
        hess_expected[1] = hessian.sparseView();
    }

    SECTION("Vertical Displ Large Gap")
    {
        vertices.row(2) << 0.0, 0.5;
        vertices.row(3) << 0.0, 1.0;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;

        hess_expected.resize(2);
        hessian << 1.8844221694086e-15, -1.8844221694086e-15, 0, 0,
            19.9070341505171, -19.9070341505171, 3.7688443388172e-15, 0,
            -1.8844221694086e-15, 1.8844221694086e-15, 0, 0, 19.9070341505171,
            -19.9070341505171, -3.7688443388172e-15, 0, 0, 0, 0, 0,
            -39.8140683010341, 39.8140683010341, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            19.9070341505171, 19.9070341505171, -39.8140683010341, 0,
            373.01213170559, 373.01213170559, -746.02426341118, 0,
            -19.9070341505171, -19.9070341505171, 39.8140683010341, 0,
            373.01213170559, 373.01213170559, -746.02426341118, 0,
            3.7688443388172e-15, -3.7688443388172e-15, 0, 0, -746.02426341118,
            -746.02426341118, 1492.04852682236, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hess_expected.resize(2);
        hessian << 9.42211084704301e-16, -9.42211084704301e-16, 0, 0,
            9.95351707525853, -9.95351707525853, 1.8844221694086e-15, 0,
            -9.42211084704301e-16, 9.42211084704301e-16, 0, 0, 9.95351707525853,
            -9.95351707525853, -1.8844221694086e-15, 0, 0, 0, 0, 0,
            -19.9070341505171, 19.9070341505171, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            9.95351707525853, 9.95351707525853, -19.9070341505171, 0,
            186.506065852795, 186.506065852795, -373.01213170559, 0,
            -9.95351707525853, -9.95351707525853, 19.9070341505171, 0,
            186.506065852795, 186.506065852795, -373.01213170559, 0,
            1.8844221694086e-15, -1.8844221694086e-15, 0, 0, -373.01213170559,
            -373.01213170559, 746.02426341118, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[1] = hessian.sparseView();
    }

    SECTION("Horizontal Displ Small Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.499;
        displacements.row(3) << 0.0, -0.499;

        hess_expected.resize(4);
        hessian << 0, 0, 0, 0, 66933.4186463403, -66933.4186463403, 0, 0, 0, 0,
            0, 0, 16733.3546615851, -16733.3546615851, 0, 0, 0, 0, 0, 0,
            -83666.7733079253, 83666.7733079253, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            66933.4186463403, 16733.3546615851, -83666.7733079253, 0,
            106880432.109393, 26720108.0273483, -133600540.136742, 0,
            -66933.4186463403, -16733.3546615851, 83666.7733079253, 0,
            26720108.0273483, 6680027.00683709, -33400135.0341854, 0, 0, 0, 0,
            0, -133600540.136742, -33400135.0341854, 167000675.170927, 0, 0, 0,
            0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 0, 0, 0, 0, 40160.0511878042, -40160.0511878042, 0, 0, 0, 0,
            0, 0, 10040.012796951, -10040.012796951, 0, 0, 0, 0, 0, 0,
            -50200.0639847552, 50200.0639847552, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            40160.0511878042, 10040.012796951, -50200.0639847552, 0,
            64128259.265636, 16032064.816409, -80160324.082045, 0,
            -40160.0511878042, -10040.012796951, 50200.0639847552, 0,
            16032064.816409, 4008016.20410225, -20040081.0205113, 0, 0, 0, 0, 0,
            -80160324.082045, -20040081.0205113, 100200405.102556, 0, 0, 0, 0,
            0, 0, 0, 0, 0;
        hess_expected[1] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 0, 0, 0, 0, 16733.3546615851, -16733.3546615851, 0, 0, 0, 0,
            0, 0, 66933.4186463403, -66933.4186463403, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, -83666.7733079253, 83666.7733079253, 0, 0,
            16733.3546615851, 66933.4186463403, 0, -83666.7733079253,
            6680027.00683709, 26720108.0273483, 0, -33400135.0341854,
            -16733.3546615851, -66933.4186463403, 0, 83666.7733079253,
            26720108.0273483, 106880432.109393, 0, -133600540.136742, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -33400135.0341854, -133600540.136742, 0,
            167000675.170927;
        hess_expected[2] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 0, 0, 0, 0, 10040.012796951, -10040.012796951, 0, 0, 0, 0, 0,
            0, 40160.0511878042, -40160.0511878042, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, -50200.0639847552, 50200.0639847552, 0, 0,
            10040.012796951, 40160.0511878042, 0, -50200.0639847552,
            4008016.20410225, 16032064.816409, 0, -20040081.0205113,
            -10040.012796951, -40160.0511878042, 0, 50200.0639847552,
            16032064.816409, 64128259.265636, 0, -80160324.082045, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, -20040081.0205113, -80160324.082045, 0,
            100200405.102556;
        hess_expected[3] = hessian.sparseView();
    }

    SECTION("Horizontal Displ Big Gap")
    {
        vertices.row(2) << -0.3, 0.5;
        vertices.row(3) << 0.3, 0.5;
        displacements.row(2) << 0.0, -0.45;
        displacements.row(3) << 0.0, -0.45;

        hess_expected.resize(4);
        hessian << 1.8844221694086e-15, -1.8844221694086e-15, 0, 0,
            31.8512546408273, -31.8512546408273, 3.7688443388172e-15, 0,
            -1.8844221694086e-15, 1.8844221694086e-15, 0, 0, 7.96281366020682,
            -7.96281366020683, -3.7688443388172e-15, 0, 0, 0, 0, 0,
            -39.8140683010341, 39.8140683010341, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            31.8512546408273, 7.96281366020682, -39.8140683010341, 0,
            954.91105716631, 238.727764291578, -1193.63882145789, 0,
            -31.8512546408273, -7.96281366020683, 39.8140683010341, 0,
            238.727764291578, 59.6819410728944, -298.409705364472, 0,
            3.7688443388172e-15, -3.7688443388172e-15, 0, 0, -1193.63882145789,
            -298.409705364472, 1492.04852682236, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        hess_expected[0] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 1.13065330164516e-15, -1.13065330164516e-15, 0, 0,
            19.1107527844964, -19.1107527844964, 2.26130660329032e-15, 0,
            -1.13065330164516e-15, 1.13065330164516e-15, 0, 0, 4.7776881961241,
            -4.7776881961241, -2.26130660329032e-15, 0, 0, 0, 0, 0,
            -23.8884409806205, 23.8884409806205, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            19.1107527844964, 4.7776881961241, -23.8884409806205, 0,
            572.946634299786, 143.236658574947, -716.183292874733, 0,
            -19.1107527844964, -4.7776881961241, 23.8884409806205, 0,
            143.236658574947, 35.8091646437367, -179.045823218683, 0,
            2.26130660329032e-15, -2.26130660329032e-15, 0, 0,
            -716.183292874733, -179.045823218683, 895.229116093416, 0, 0, 0, 0,
            0, 0, 0, 0, 0;
        hess_expected[1] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 1.8844221694086e-15, -1.8844221694086e-15, 0, 0,
            7.96281366020683, -7.96281366020682, 0, 3.7688443388172e-15,
            -1.8844221694086e-15, 1.8844221694086e-15, 0, 0, 31.8512546408273,
            -31.8512546408273, 0, -3.7688443388172e-15, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, -39.8140683010341, 39.8140683010341, 0, 0,
            7.96281366020683, 31.8512546408273, 0, -39.8140683010341,
            59.6819410728944, 238.727764291578, 0, -298.409705364472,
            -7.96281366020682, -31.8512546408273, 0, 39.8140683010341,
            238.727764291578, 954.91105716631, 0, -1193.63882145789, 0, 0, 0, 0,
            0, 0, 0, 0, 3.7688443388172e-15, -3.7688443388172e-15, 0, 0,
            -298.409705364472, -1193.63882145789, 0, 1492.04852682236;
        hess_expected[2] = hessian.sparseView();
        hess_expected.resize(4);
        hessian << 1.13065330164516e-15, -1.13065330164516e-15, 0, 0,
            4.7776881961241, -4.7776881961241, 0, 2.26130660329032e-15,
            -1.13065330164516e-15, 1.13065330164516e-15, 0, 0, 19.1107527844964,
            -19.1107527844964, 0, -2.26130660329032e-15, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, -23.8884409806205, 23.8884409806205, 0, 0,
            4.7776881961241, 19.1107527844964, 0, -23.8884409806205,
            35.8091646437367, 143.236658574947, 0, -179.045823218683,
            -4.7776881961241, -19.1107527844964, 0, 23.8884409806205,
            143.236658574947, 572.946634299786, 0, -716.183292874733, 0, 0, 0,
            0, 0, 0, 0, 0, 2.26130660329032e-15, -2.26130660329032e-15, 0, 0,
            -179.045823218683, -716.183292874733, 0, 895.229116093416;
        hess_expected[3] = hessian.sparseView();
    }

    REQUIRE(
        barrier.initial_epsilon == ccd::opt::InitialBarrierEpsilon::MIN_TOI);
    barrier.initialize(vertices, edges, displacements * 2.0);

    barrier.compute_constraints_hessian(displacements, hess_actual);
    CHECK(hess_actual.size() == hess_expected.size());

    for (size_t i = 0; i < hess_actual.size(); i++) {
        CHECK(
            (hess_actual[i] - hess_expected[i]).toDense().squaredNorm() < 1e-6);
    }

    Eigen::VectorXd v_actual;
    Eigen::MatrixXd jac_actual;

    barrier.compute_constraints_and_derivatives(
        displacements, v_actual, jac_actual, hess_actual);
    for (size_t i = 0; i < hess_actual.size(); i++) {
        CHECK(
            (hess_actual[i] - hess_expected[i]).toDense().squaredNorm() < 1e-6);
    }
}
