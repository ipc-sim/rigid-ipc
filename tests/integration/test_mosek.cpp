#ifdef BUILD_WITH_MOSEK

#include <Eigen/SparseCore>
#include <catch2/catch.hpp>

#include <igl/mosek/mosek_quadprog.h>

#define INF_D (std::numeric_limits<double>::infinity())

TEST_CASE("Simple tests of MOSEK", "[opt][mosek]")
{
    // Load problem data
    const int N = 4;
    const int M = N;

    // Quadratic matrix
    Eigen::SparseMatrix<double> Q(N, N);
    Q.setIdentity();

    // Quadratic linear term
    Eigen::VectorXd c(N);
    c.setZero();

    // Quadratic constant term
    double cf = 0.0;

    // Linear constraint matrix
    Eigen::SparseMatrix<double> A(M, N);
    A.setIdentity();

    // Linear constraint lower bounds
    Eigen::VectorXd lc(M);
    lc.setOnes();
    Eigen::VectorXd expected_solution(M);
    expected_solution.setOnes();
    SECTION("No Lower Bounds")
    {
        lc *= -INF_D;
        expected_solution *= 0;
    }
    SECTION("Lower Bounds = 1") {}
    SECTION("Lower Bounds = 10")
    {
        lc *= 10;
        expected_solution *= 10;
    }
    // Linear constraint upper bounds
    Eigen::VectorXd uc(M);
    uc.setOnes();
    uc *= INF_D;
    Eigen::VectorXd lx(M);
    lx.setOnes();
    lx *= -INF_D;
    Eigen::VectorXd ux(M);
    ux.setOnes();
    ux *= INF_D;

    igl::mosek::MosekData mosek_data;

    Eigen::VectorXd x(N);
    bool success = igl::mosek::mosek_quadprog(
        Q, c, cf, A, lc, uc, lx, ux, mosek_data, x);

    // I always require success in everything I do!
    REQUIRE(success);
    REQUIRE(x.size() == expected_solution.size());
    CHECK((x - expected_solution).squaredNorm() == Approx(0.0).margin(1e-12));
}

#endif
