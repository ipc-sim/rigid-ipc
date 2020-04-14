#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <autodiff/autodiff_types.hpp>
#include <opt/constrained_problem.hpp>
#include <solvers/ncp_solver.hpp>

#include <logger.hpp>

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

    SECTION("Abs Value Case")
    {

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

    SECTION("Circle Case")
    {
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

    class AdHocProblem : public virtual ConstrainedProblem {
    public:
        AdHocProblem(
            Eigen::SparseMatrix<double>& _A,
            Eigen::VectorXd& _b,
            std::function<DVector(const Eigen::VectorXd& x)>& _gdiff)
            : A(_A)
            , b(_b)
            , gdiff(_gdiff)
        {
            is_dof_fixed_ = Eigen::VectorXb::Zero(NUM_VARS);
        }

        virtual void compute_objective(
            const Eigen::VectorXd& x,
            double& fx,
            Eigen::VectorXd& grad_fx,
            Eigen::SparseMatrix<double>& hess_fx,
            bool compute_grad = true,
            bool compute_hess = true) override
        {
            grad_fx = A * x - b;
            fx = (grad_fx).squaredNorm() / 2.0;
            if (compute_hess) {
                hess_fx = A;
            }
        }

        virtual void compute_constraints(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& jac_gx,
            std::vector<Eigen::SparseMatrix<double>>& hess_gx,
            bool compute_grad = true,
            bool compute_hess = true) override
        {
            DVector g = gdiff(x);
            gx = Eigen::VectorXd(g.rows());
            jac_gx = Eigen::MatrixXd(gx.rows(), NUM_VARS);
            for (int i = 0; i < g.rows(); ++i) {
                gx(i) = g(i).getValue();
                jac_gx.row(i) = g(i).getGradient();
            }

            if (compute_hess) {
                throw "not computing hess_gx";
            }
        }

        using ConstrainedProblem::compute_constraints;

        void compute_constraints_using_normals(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& jac_gx) override
        {
            return compute_constraints(x, gx, jac_gx);
        }

        const Eigen::VectorXd& starting_point() override { return b; }
        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }

        virtual int num_vars() const override { return NUM_VARS; }

        virtual bool has_collisions(
            const Eigen::VectorXd& xi, const Eigen::VectorXd& xj) const override
        {
            return false;
        }
        virtual double
        compute_min_distance(const Eigen::VectorXd& x) const override
        {
            return -1;
        };

        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        std::function<DVector(const Eigen::VectorXd& x)> gdiff;
        Eigen::VectorXb is_dof_fixed_;
    };

    Eigen::VectorXd x(NUM_VARS), alpha(NUM_CONSTRAINTS);
    AdHocProblem problem(A, b, g_diff);

    NCPSolver solver;
    solver.set_problem(problem);
    solver.max_iterations = 3000;
    solver.convergence_tolerance = 1E-8;
    solver.do_line_search = false;
    solver.solve_for_active_cstr = false;
    solver.update_type = NCPUpdate::G_GRADIENT;

    OptimizationResults results;

    // Solve using Guass-Seidel
    solver.lcp_solver = LCPSolver::LCP_GAUSS_SEIDEL;
    results = solver.solve();
    CHECK(results.finished);
    CHECK(results.success);
    CHECK((expected - results.x).squaredNorm() < 1E-6);

    // Solve using Fischer-Newton
    solver.lcp_solver = LCPSolver::LCP_NEWTON;
    results = solver.solve();
    CHECK(results.finished);
    CHECK(results.success);
    CHECK((expected - results.x).squaredNorm() < 1E-6);

#ifdef BUILD_WITH_MOSEK
    // Solve using Mosek QP
    solver.lcp_solver = LCPSolver::LCP_MOSEK;
    results = solver.solve();
    CHECK(results.finished);
    CHECK(results.success);
    CHECK((expected - results.x).squaredNorm() < 1E-6);
#endif
}
