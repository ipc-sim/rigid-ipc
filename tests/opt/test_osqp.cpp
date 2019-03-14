#ifdef BUILD_WITH_OSQP

#include <Eigen/Sparse>
#include <catch.hpp>
#include <osqp.h>

#define INF_D (std::numeric_limits<double>::infinity())

TEST_CASE("Simple tests of OSQP", "[opt][nlopt]")
{
    // Load problem data
    const c_int N = 4;
    const c_int M = N;

    // Quadratic matrix
    Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> P(N, N);
    P.setIdentity();

    // Quadratic linear term
    Eigen::Matrix<c_float, N, 1> q;
    q.setZero();

    // Linear constraint matrix
    Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> A(M, N);
    A.setIdentity();

    // Linear constraint lower bounds
    Eigen::Matrix<c_float, M, 1> l;
    l.setOnes();
    Eigen::Matrix<c_float, M, 1> expected_solution;
    expected_solution.setOnes();
    SECTION("No Lower Bounds")
    {
        l *= -INF_D;
        expected_solution *= 0;
    }
    SECTION("Lower Bounds = 1") {}
    SECTION("Lower Bounds = 10")
    {
        l *= 10;
        expected_solution *= 10;
    }
    // Linear constraint upper bounds
    Eigen::Matrix<c_float, M, 1> u;
    u.setOnes();
    u *= INF_D;

    // Populate data
    OSQPData data; // OSQPData
    data.n = N;
    data.m = M;
    data.P = csc_matrix(P.rows(), P.cols(), P.nonZeros(), P.valuePtr(),
        P.innerIndexPtr(), P.outerIndexPtr());
    data.q = q.data();
    data.A = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(),
        A.innerIndexPtr(), A.outerIndexPtr());
    data.l = l.data();
    data.u = u.data();

    // Define Solver settings as default
    OSQPSettings settings;
    osqp_set_default_settings(&settings);
    settings.alpha = 1.0; // Change alpha parameter

    // Setup workspace
    OSQPWorkspace* work(osqp_setup(&data, &settings)); // Workspace

    // Solve Problem
    osqp_solve(work);

    Eigen::Map<Eigen::Matrix<c_float, N, 1>> x(work->solution->x, N);

    REQUIRE(x.size() == expected_solution.size());
    CHECK((x - expected_solution).squaredNorm() == Approx(0.0).margin(1e-12));

    // Cleanup
    osqp_cleanup(work);
}

#endif
