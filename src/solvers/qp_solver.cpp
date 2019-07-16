#include <solvers/qp_solver.hpp>

#include <solvers/optimization_solver.hpp>

#include <Eigen/SparseCore>
#ifdef BUILD_WITH_OSQP
#include <osqp.h>
#endif
#if BUILD_WITH_MOSEK
#include <igl/mosek/mosek_quadprog.h>
#endif

#include <utils/not_implemented_error.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    QPSolver::QPSolver()
        : absolute_tolerance(1e-8)
        , relative_tolerance(1e-8)
        , qp_solver(QPSolverType::OSQP)
    {
    }
    QPSolver::~QPSolver() {}

    void QPSolver::quadratic_energy(OptimizationProblem& problem)
    {
        // Approximate to Quadratic Energy of the form
        //  ½ * x^T * Q * x + c^T * x + cf

        Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.num_vars);
        Q = problem.eval_hessian_f(x0);
        c = problem.eval_grad_f(x0);
        cf = problem.eval_f(x0);
    }

    void QPSolver::linearized_constraints(OptimizationProblem& problem)
    {
        // Linearized constraints
        // 0.0 ≤ g(x0) + ∇g(x0)(x - x0) ≤ INF →
        // 0.0 - g(x0) + ∇g(x0) * x0 ≤ ∇g(x0) * x ≤ INF
        // Linear constraint matrix
        // A = ∇g(x0) ∈ R^(m × n)
        A = problem.eval_jac_g(problem.x0).sparseView();
        // Linear constraint lower bounds
        // (0.0 - g(x0) + ∇g(x0) * x0) ∈ R^m
        lc = - problem.eval_g(problem.x0) + A * problem.x0;
        // Linear constraint upper bounds
        // u ∈ R^m
        uc = Eigen::VectorXd(lc.rows());
        uc.setConstant(NO_UPPER_BOUND);
    }

    OptimizationResults QPSolver::solve(OptimizationProblem& problem)
    {
        quadratic_energy(problem);
        linearized_constraints(problem);

        OptimizationResults results;
        results.x.resize(problem.num_vars, 1);

        switch (qp_solver) {
        case QPSolverType::OSQP:
            results.success = solve_with_osqp(results.x);
            break;
        case QPSolverType::MOSEK:
            results.success = solve_with_mosek(results.x);
            break;
        }

        // Compute objective at x
        results.minf = problem.eval_f(results.x);

        // Check the solve was successful
        auto lhs = (A * results.x).array();
        results.success &= results.minf >= 0;

        return results;
    }

    bool QPSolver::solve_with_mosek(Eigen::MatrixXd& x)
    {
#ifdef BUILD_WITH_MOSEK
        igl::mosek::MosekData mosek_data;
        Eigen::VectorXd y;
        bool success = igl::mosek::mosek_quadprog(
            Q, c, cf, A, lc, uc, lx, ux, mosek_data, y);
        x << y;
        return success;

#else
        throw NotImplementedError("MOSEK is not enabled");
#endif
    }

    bool QPSolver::solve_with_osqp(Eigen::MatrixXd& x)
    {

#ifdef BUILD_WITH_OSQP
        int n = int(c.rows()), m = int(lc.rows());
        assert(Q.rows() == n && Q.cols() == n);
        assert(A.rows() == m && A.cols() == n);
        assert(uc.rows() == m);
        assert(x.size() == n);

        // Populate data
        OSQPData data; // OSQPData
        data.n = n;
        data.m = m;
        Q.makeCompressed();
        data.P = csc_matrix(Q.rows(), Q.cols(), Q.nonZeros(), Q.valuePtr(),
            Q.innerIndexPtr(), Q.outerIndexPtr());
        data.q = c.data();
        A.makeCompressed();
        data.A = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(),
            A.innerIndexPtr(), A.outerIndexPtr());
        data.l = lc.data();
        data.u = uc.data();

        // Define Solver settings as default
        OSQPSettings osqp_settings;
        osqp_set_default_settings(&osqp_settings);
        osqp_settings.max_iter = max_iterations;
        osqp_settings.eps_abs = absolute_tolerance;
        osqp_settings.eps_rel = relative_tolerance;
        osqp_settings.eps_prim_inf = relative_tolerance / 10;
        osqp_settings.eps_dual_inf = relative_tolerance / 10;
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
#else
        throw NotImplementedError("OSQP is not enabled");
#endif
    }

} // namespace opt
} // namespace ccd
