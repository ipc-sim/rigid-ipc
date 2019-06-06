#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <autodiff/autodiff_types.hpp>
#include <opt/ncp_solver.hpp>

// ---------------------------------------------------
// SETUP
// ---------------------------------------------------
static const int NUM_VARS = 2;
static const int NUM_CONSTRAINTS = 2;

// differentiable helpers
template <typename T> using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;

typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef DScalar1<double, Eigen::Matrix<double, NUM_VARS, 1>> DScalar;
typedef VectorXT<DScalar> DVector;

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("NCP", "[opt][NCP][NCP-Interface]")
{
    using namespace ccd::opt;
    DiffScalarBase::setVariableCount(size_t(NUM_VARS));

    Eigen::SparseMatrix<double> A(NUM_VARS, NUM_VARS);
    Eigen::VectorXd b(NUM_VARS), expected(NUM_VARS);

    std::function<DVector(const Eigen::VectorXd& x)> g_diff;

    // ------------------------------------------------------------------------
    // PROBLEM SETUP
    // ------------------------------------------------------------------------
    A.setIdentity();
    b << -1, -2.5;

    SECTION("Linear Case")
    {
        g_diff = [](const Eigen::VectorXd& x) -> DVector {
            DVector gx(NUM_CONSTRAINTS);
            DScalar x0(0, x[0]);
            DScalar x1(1, x[1]);

            gx(0) = x0;
            gx(1) = x1;
            return gx;
        };
        expected << 0.0, 0.0;
    }

    SECTION("Quadratic Case")
    {
        g_diff = [](const Eigen::VectorXd& x) -> DVector {
            DVector gx(NUM_CONSTRAINTS);
            DScalar x0(0, x[0]);
            DScalar x1(1, x[1]);

            gx(0) = 0.04 - x0 * x0;
            gx(1) = 0.09 - x1 * x1;
            return gx;
        };
        expected << -0.2, -0.3;
    }

    SECTION("Abs Value Case"){

        g_diff = [](const Eigen::VectorXd& x) -> DVector {
            DVector gx(NUM_CONSTRAINTS);
            DScalar x0(0, x[0]);
            DScalar x1(1, x[1]);

            gx(0) = 0.2 - (x0 > 0 ? x0 : -x0);
            gx(1) = 0.3 - (x1 > 0 ? x1 : -x1);
            return gx;
        };
        expected << -0.2, -0.3;
    }

    SECTION("Ciecle Case"){

        g_diff = [](const Eigen::VectorXd& x) -> DVector {
            DVector gx(NUM_CONSTRAINTS);
            DScalar x0(0, x[0]);
            DScalar x1(1, x[1]);

            gx(0) = 1.0 - (x0 - 1.0) * (x0 - 1.0);
            gx(1) = 1.0 - (x1 - 2.5) * (x1 - 2.5);
            return gx;
        };
        expected << 0.0, 1.5;
    }
    AdHocProblem problem;
    problem.g = [&g_diff](const Eigen::VectorXd x) -> Eigen::VectorXd {
        DVector gx = g_diff(x);
        Eigen::VectorXd g(gx.rows());
        for (int i = 0; i < gx.rows(); ++i) {
            g(i) = gx(i).getValue();
        }

        return g;
    };

    problem.jac_g = [&g_diff](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        DVector gx = g_diff(x);
        Eigen::MatrixXd jac_gx(gx.rows(), NUM_VARS);
        for (int i = 0; i < gx.rows(); ++i) {
            jac_gx.row(i) = gx(i).getGradient();
        }

        return jac_gx;
    };

    Eigen::VectorXd x(NUM_VARS), alpha(NUM_CONSTRAINTS);
    NCPSolver solver;
    solver.max_iterations= 300;
    solver.convergence_tolerance  = 1E-8;
    solver.keep_in_unfeasible = false;
    solver.check_convergence = false;
    solver.solve_for_active_cstr = false;
    solver.update_type = NcpUpdate::LINEARIZED;
    solver.lcp_solver = LCPSolver::LCP_GAUSS_SEIDEL;

    bool success = solver.solve_ncp(A, b, problem, x, alpha);

    CHECK(success);
    CHECK((expected - x).squaredNorm() < 1E-6);
}
