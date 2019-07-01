#ifdef BUILD_WITH_OSQP

#include <Eigen/Sparse>
#include <catch2/catch.hpp>
#include <iostream>
#include <osqp.h>

TEST_CASE("OSQP_integration", "[opt][osqp][!mayfail]")
{
    typedef Eigen::Matrix<c_float, Eigen::Dynamic, 1> VectorX_OSQP;
    typedef Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int>
        SparseMatrix_OSQP;

    // Load problem data
    const c_int N = 4;
    const c_int M = N;

    // Quadratic matrix
    SparseMatrix_OSQP P(N, N);
    P.setIdentity();

    // Quadratic linear term
    VectorX_OSQP q = VectorX_OSQP::Zero(N);

    // Linear constraint matrix
    SparseMatrix_OSQP A(M, N);
    A.setIdentity();

    // Linear constraint lower bounds
    VectorX_OSQP l(M);
    VectorX_OSQP expected_solution(M);
    SECTION("No Lower Bounds")
    {
        l.setConstant(-2e19);
        expected_solution.setZero();
    }
    SECTION("Lower Bounds = 1")
    {
        l.setOnes();
        expected_solution.setOnes();
    }
    SECTION("Lower Bounds = 10")
    {
        l.setConstant(10);
        expected_solution.setConstant(10);
    }
    // Linear constraint upper bounds
    VectorX_OSQP u = VectorX_OSQP::Constant(M, 2e19);

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
    OSQPWorkspace* work = osqp_setup(&data, &settings); // Workspace

    // Solve Problem
    osqp_solve(work);

    Eigen::Map<VectorX_OSQP> x(work->solution->x, N);

    REQUIRE(x.size() == expected_solution.size());
    CHECK((x - expected_solution).squaredNorm() == Approx(0.0).margin(1e-8));
    // Cleanup
    osqp_cleanup(work);
}

#endif
