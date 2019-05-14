#include <opt/linearized_constraint_solver.hpp>
#include <opt/solver.hpp>

#include <Eigen/SparseCore>
#ifdef BUILD_WITH_OSQP
#include <osqp.h>
#endif
#if BUILD_WITH_MOSEK
#include <igl/mosek/mosek_quadprog.h>
#endif

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>

#include <iostream>

namespace ccd {
namespace opt {

#ifdef BUILD_WITH_OSQP
    bool solve_qp_with_osqp(
        Eigen::SparseMatrix<double, Eigen::ColMajor, int>& P,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& q,
        Eigen::SparseMatrix<double, Eigen::ColMajor, int>& A,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& l,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& u, Eigen::MatrixXd& x,
        const SolverSettings& settings)
    {
        int n = int(q.rows()), m = int(l.rows());
        assert(P.rows() == n && P.cols() == n);
        assert(A.rows() == m && A.cols() == n);
        assert(u.rows() == m);
        assert(x.size() == n);

        // Populate data
        OSQPData data; // OSQPData
        data.n = n;
        data.m = m;
        P.makeCompressed();
        data.P = csc_matrix(P.rows(), P.cols(), P.nonZeros(), P.valuePtr(),
            P.innerIndexPtr(), P.outerIndexPtr());
        data.q = q.data();
        A.makeCompressed();
        data.A = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(),
            A.innerIndexPtr(), A.outerIndexPtr());
        data.l = l.data();
        data.u = u.data();

        // Define Solver settings as default
        OSQPSettings osqp_settings;
        osqp_set_default_settings(&osqp_settings);
        osqp_settings.max_iter = settings.max_iter;
        osqp_settings.eps_abs = settings.absolute_tolerance;
        osqp_settings.eps_rel = settings.relative_tolerance;
        osqp_settings.eps_prim_inf = settings.relative_tolerance / 10;
        osqp_settings.eps_dual_inf = settings.relative_tolerance / 10;
        osqp_settings.polish = true;

        // Setup workspace
        OSQPWorkspace* work(osqp_setup(&data, &osqp_settings)); // Workspace

        // Solve Problem
        int exit_flag = osqp_solve(work);

        // Save solution
        x = Eigen::Map<Eigen::VectorXd>(work->solution->x, n);

        // Cleanup
        osqp_cleanup(work);
        return exit_flag >= 0;
    }
#endif

    // Optimize the displacments using linearized constraints
    OptimizationResults solve_problem_with_linearized_constraints(
        OptimizationProblem& problem, const SolverSettings& settings)
    {
        // Quadratic Energy
        // 1/2 * || U - U0 ||^2 = (U - U0)^T * (U - U0) =
        //      1/2 * U^T * I * U - U0^T * U + 1/2 * U0^T * U0

        // Quadratic matrix P = I
        Eigen::SparseMatrix<double> P(problem.num_vars, problem.num_vars);
        P.setIdentity();

        // Linear term
        // q = -U0
        Eigen::VectorXd q = -problem.x0;

        // Constant term
        // c = 1/2 * U0^T * U0
        double c = 0.5 * problem.x0.transpose() * problem.x0;

        // Check that the QP is properly expressed
        assert((((0.5 * problem.x0.transpose() * P * problem.x0
                     + q.transpose() * problem.x0)
                        .array()
                    + c)
                    .abs()
            < 1e-8)
                   .all());

        // Linearized constraints
        // ℓ ≤ g(x0) + ∇g(x0)(x - x0) ≤ u →
        // ℓ - g(x0) + ∇g(x0) * x0 ≤ ∇g(x0) * x ≤ u
        // Linear constraint matrix
        // A = ∇g(x0) ∈ R^(m × n)
        Eigen::SparseMatrix<double> A = problem.eval_jac_g(problem.x0).sparseView();
        // Linear constraint lower bounds
        // (ℓ - g(x0) + ∇g(x0) * x0) ∈ R^m
        Eigen::VectorXd l
            = problem.g_lower - problem.eval_g(problem.x0) + A * problem.x0;
        // Linear constraint upper bounds
        // u ∈ R^m
        Eigen::VectorXd u = problem.g_upper;

        OptimizationResults results;
        results.x.resize(problem.num_vars, 1);

        switch (settings.qp_solver) {
        case OSQP:
#ifdef BUILD_WITH_OSQP
            results.success
                = solve_qp_with_osqp(P, q, A, l, u, results.x, settings);
            break;
#else
            throw NotImplementedError("OSQP is not enabled");
#endif
        case MOSEK:
#ifdef BUILD_WITH_MOSEK
            igl::mosek::MosekData mosek_data;
            Eigen::VectorXd x(problem.num_vars);
            results.success = igl::mosek::mosek_quadprog(P, q, c, A, l, u,
                problem.x_lower, problem.x_upper, mosek_data, x);
            results.x = x;
            break;
#else
            throw NotImplementedError("MOSEK is not enabled");
#endif
        }

        // Compute objective at x
        results.minf = problem.eval_f(results.x);

        // Check the solve was successful
        auto lhs = (A * results.x).array();
        results.success &= results.minf >= 0
            && ((l.array() - 10 * settings.relative_tolerance) <= lhs).all()
            && (lhs <= (u.array() + 10 * settings.relative_tolerance)).all();

        return results;
    }

} // namespace opt
} // namespace ccd
