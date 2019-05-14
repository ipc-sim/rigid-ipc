// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/displacement_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>

#include <autodiff/finitediff.hpp>

#include <logger.hpp>

#include <profiler.hpp>

namespace ccd {
namespace opt {

    // Optimize the displacment opt problem with the given method and starting
    // value.
    OptimizationResults displacement_optimization(OptimizationProblem& problem,
        const Eigen::MatrixX2d& U0, SolverSettings& settings)
    {
#ifdef PROFILE_FUNCTIONS
        reset_profiler();
        igl::Timer timer;
        timer.start();
#endif

        // initial value
        Eigen::MatrixXd x0 = U0;
        x0.resize(U0.size(), 1); // Flatten displacements
        problem.x0 = x0;

        OptimizationResults result;
        result = solve_problem(problem, settings);
        result.x.resize(U0.rows(), 2); // Unflatten displacments

#ifdef PROFILE_FUNCTIONS
        timer.stop();
        print_profile(timer.getElapsedTime());
#endif

        return result;
    } // namespace opt

    ParticlesDisplProblem::ParticlesDisplProblem()
        : constraint(nullptr)
    {
    }

    ParticlesDisplProblem::~ParticlesDisplProblem() {}

    void ParticlesDisplProblem::initialize(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2i& E, const Eigen::MatrixX2d& U,
        CollisionConstraint& cstr)
    {
        vertices = V;
        edges = E;
        displacements = U;

        constraint = &cstr;
        constraint->initialize(V, E, U);
        initProblem();
    }

    void ParticlesDisplProblem::initProblem()
    {
        u_ = displacements;
        u_.resize(displacements.size(), 1);

        num_vars = int(u_.size());
        x0.resize(num_vars);
        x_lower.resize(num_vars);
        x_upper.resize(num_vars);
        x_lower.setConstant(NO_LOWER_BOUND);
        x_upper.setConstant(NO_UPPER_BOUND);

        // TODO: this is wrong, it will depend on constraint type
        num_constraints = int(edges.rows());
        g_lower.resize(num_constraints);
        g_upper.resize(num_constraints);
        g_lower.setConstant(0.0);
        g_upper.setConstant(NO_UPPER_BOUND);
    }

    double ParticlesDisplProblem::eval_f(const Eigen::VectorXd& x)
    {
        return (x - u_).squaredNorm() / 2.0;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return (x - u_);
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        return Eigen::MatrixXd::Identity(x.size(), x.size());
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        Eigen::VectorXd gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints(vertices, edges, Uk, gx);
        return gx;
    };

    Eigen::MatrixXd ParticlesDisplProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        Eigen::MatrixXd jac_gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints_jacobian(vertices, edges, Uk, jac_gx);

        return jac_gx;
    };

    std::vector<Eigen::MatrixXd> ParticlesDisplProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        std::vector<Eigen::MatrixXd> hess_gx;
        if (constraint->recompute_collision_set) {
            constraint->detecteCollisions(vertices, edges, Uk);
        }
        constraint->compute_constraints_hessian(vertices, edges, Uk, hess_gx);
        return hess_gx;
    };

    NCPDisplacementOptimization::NCPDisplacementOptimization()
        : max_iterations(100)
        , update_method(NcpUpdate::LINEARIZED)
        , lcp_solver(LCPSolver::LCP_GAUSS_SEIDEL)
        , keep_in_unfeasible(true)
        , check_convergence(false)
        , convegence_tolerance(1e-6)
    {
    }

    OptimizationResults NCPDisplacementOptimization::solve(OptimizationProblem& problem)
    {
        // Solves the KKT conditions of the Optimization Problem
        //  (U - Uk) = \nabla g(U)
        //  s.t V(U) >= 0
        int num_vars = int(problem.x0.rows());
        Eigen::SparseMatrix<double> A(num_vars, num_vars);
        A.setIdentity();

        Eigen::VectorXd b, x_opt, lambda_opt;
        // obtain Uk from grad_f using x=0
        b = -problem.eval_grad_f(Eigen::VectorXd::Zero(num_vars));

        int num_it = 0;
        callback_intermediate_ncp callback
            = [&](const Eigen::VectorXd& x, const Eigen::VectorXd& alpha,
                  const double gamma) {
                  // settings->intermediate_cb(x, problem->f(x), alpha, gamma,
                  // num_it);
                  num_it += 1;
              };
        OptimizationResults result;
        result.success = solve_ncp(A, b, problem, max_iterations, callback,
            update_method, lcp_solver, x_opt, lambda_opt, keep_in_unfeasible,
            check_convergence, convegence_tolerance);
        result.x = x_opt;
        result.minf = problem.eval_f(x_opt);

        return result;
    }

} // namespace opt
} // namespace ccd
