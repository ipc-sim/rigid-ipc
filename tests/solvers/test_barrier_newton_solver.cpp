#include <catch2/catch.hpp>
#include <cmath>

#include <finitediff.hpp>

#include <autodiff/autodiff_types.hpp>
#include <barrier/barrier.hpp>
#include <solvers/homotopy_solver.hpp>
#include <solvers/newton_solver.hpp>
#include <utils/not_implemented_error.hpp>

#include <logger.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE(
    "Check barrier problem derivatives", "[opt][barrier][barrier_problem]")
{
    // TODO: Generate a random problem and test the derivatives
    int num_vars = GENERATE(1, 2, 5, 10);

    class AdHocProblem : public virtual BarrierProblem {
    public:
        AdHocProblem(int num_vars, double epsilon)
            : num_vars_(num_vars)
            , m_barrier_epsilon(epsilon)
        {
        }

        double compute_energy_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad_Ex,
            Eigen::SparseMatrix<double>& hess_Ex,
            bool compute_grad = true,
            bool compute_hess = true) override
        {
            if (compute_grad) {
                grad_Ex = x;
            }
            if (compute_hess) {
                hess_Ex =
                    Eigen::MatrixXd::Identity(x.rows(), x.rows()).sparseView();
            }
            return x.squaredNorm() / 2.0;
        }

        double compute_barrier_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad_Bx,
            Eigen::SparseMatrix<double>& hess_Bx,
            int& num_constraints,
            bool compute_grad = true,
            bool compute_hess = true) override
        {
            num_constraints = 1;

            // 1/2||x||² ≥ 1 → 1/2||x||² - 1 ≥ 0
            typedef AutodiffType<Eigen::Dynamic> Diff;
            Diff::activate(x.size());

            Eigen::Matrix<Diff::DDouble2, Eigen::Dynamic, 1> dx;

            Diff::DDouble2 b =
                poly_log_barrier(dx.squaredNorm() / 2.0 - 1, m_barrier_epsilon);

            grad_Bx = b.getGradient();
            hess_Bx = b.getHessian().sparseView();

            return b.getValue();
        }

        double barrier_hessian(double x) const override
        {
            return poly_log_barrier_hessian(x, m_barrier_epsilon);
        }

        double barrier_activation_distance() const override
        {
            return m_barrier_epsilon;
        }
        void barrier_activation_distance(const double eps) override
        {
            m_barrier_epsilon = eps;
        }

        double barrier_stiffness() const override
        {
            return m_barrier_stiffness;
        }
        void barrier_stiffness(const double kappa) override
        {
            m_barrier_stiffness = kappa;
        }

        bool
        has_collisions(const Eigen::VectorXd&, const Eigen::VectorXd&) override
        {
            return false;
        }
        virtual double compute_earliest_toi(
            const Eigen::VectorXd& xi, const Eigen::VectorXd& xj) override
        {
            return std::numeric_limits<double>::infinity();
        }
        virtual bool is_ccd_aligned_with_newton_update() override
        {
            return true;
        }

        int num_vars() const override { return num_vars_; }
        const VectorXb& is_dof_fixed() const override { return is_dof_fixed_; }
        const Eigen::VectorXd& starting_point() const { return x0; }

        double compute_min_distance(const Eigen::VectorXd& x) const override
        {
            return -1;
        }

        /// Get the world coordinates of the vertices
        Eigen::MatrixXd world_vertices(const Eigen::VectorXd& x) const override
        {
            throw NotImplementedError("no vertices");
        }

        /// Get the length of the diagonal of the worlds bounding box
        double world_bbox_diagonal() const override
        {
            throw NotImplementedError("no world bbox diagonal");
        }

        DiagonalMatrixXd mass_matrix() const override
        {
            DiagonalMatrixXd I(num_vars_);
            I.setIdentity();
            return I;
        }
        double average_mass() const override { return 1; }

        double timestep() const override { return 1; }

        int num_vars_;
        double m_barrier_epsilon;
        double m_barrier_stiffness;
        Eigen::VectorXd x0;
        VectorXb is_dof_fixed_;
    };

    double epsilon = GENERATE(1.0, 0.5, 1e-1, 5e-2);
    AdHocProblem problem(num_vars, epsilon);

    Eigen::VectorXd x(num_vars);
    x.setConstant(GENERATE(-1.0, 1.0) * GENERATE(1e-1, 1.0, 2.0 - 1e-3, 4.0));

    Eigen::VectorXd grad_fx;
    Eigen::SparseMatrix<double> hess_fx;
    double fx = problem.compute_objective(x, grad_fx, hess_fx);
    // If the function evaluates to infinity then the finite differences will
    // not work. I assume in the definition of the barrier gradient that
    // d/dx ∞ = 0.
    if (!std::isinf(fx)) {
        // Use a higher order finite difference method because the function near
        // the boundary becomes very non-linear. This problem worsens as the ϵ
        // of the boundary gets smaller.

        // Test ∇f
        Eigen::VectorXd finite_grad(problem.num_vars());
        finite_grad = ipc::rigid::eval_grad_objective_approx(problem, x);
        CHECK(fd::compare_gradient(finite_grad, grad_fx));

        // Test ∇²f
        Eigen::MatrixXd finite_hessian = eval_hess_objective_approx(problem, x);
        CHECK(fd::compare_jacobian(finite_hessian, Eigen::MatrixXd(hess_fx)));
    }
}
