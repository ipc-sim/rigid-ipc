#pragma once

#include <memory>

#include <problems/barrier_problem.hpp>
#include <solvers/newton_solver.hpp>
#include <solvers/optimization_solver.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

/**
 * @brief Solve the optimization problem using
 * Newton's Method with barriers for the constraints.
 */
class HomotopySolver : public virtual OptimizationSolver {
public:
    HomotopySolver();
    virtual ~HomotopySolver() = default;

    /// Initialize the state of the solver using the settings saved in JSON
    virtual void settings(const nlohmann::json& json) override;
    /// Export the state of the solver using the settings saved in JSON
    virtual nlohmann::json settings() const override;

    /// An identifier for the solver class
    static std::string solver_name() { return "homotopy_solver"; }
    /// An identifier for this solver
    virtual std::string name() const override
    {
        return HomotopySolver::solver_name();
    }

    virtual bool has_inner_solver() const override { return true; }

    virtual OptimizationSolver& inner_solver() override
    {
        return *inner_solver_ptr;
    }

    /// Initialize the solver with a problem to solve
    virtual void set_problem(OptimizationProblem& problem) override
    {
        assert(problem.is_barrier_problem());
        problem_ptr = dynamic_cast<BarrierProblem*>(&problem);
        inner_solver_ptr->set_problem(problem);
    }

    /// Initialize the solver state for a new solve
    virtual void init_solve(const Eigen::VectorXd& x0) override;
    /// Solve the saved optimization problem to completion
    virtual OptimizationResults solve(const Eigen::VectorXd& x0) override;
    /// Perform a single step of solving the optimization problem
    virtual OptimizationResults step_solve() override;

    ///
    /// "For a given outer iteration (i.e. fixed t) the precision of the
    /// result will be ~ m/t (in terms of how close it is to the solution of
    /// the original inequality constrained problem).  It is possible to run
    /// Newton solve to a fixed high precision (say 1e-12) to get this
    /// approximate solution, but this is wasteful.  On the other hand,
    /// running it at the same precision as the current error due to finite
    /// t is somewhat risky, especially for low values of t -- this may lead
    /// to a substantial deviation from the true solution.    So a possible
    /// approach is to use a faction c of the current outer iteration
    /// precision.  This however is excessive for high values of t,
    /// especially hitting the limits of double accuracy.  It is reasonable
    /// to limit how high we set Newton accuracy. We can try either  max(
    /// c*m/t,  e_b)  or max( c*m/t, e_max)  with e_max being the max
    /// accuracy we hope to achieve with doubles in the Newton solve (say
    /// 1e-12 or 1e-11)."
    ///
    double tinit;
    double t;
    double m;
    double c;
    double e_b;
    double t_inc;

protected:
    /// The interior newton solver of the barrier solver.
    class InnerNewtonSolver : public virtual NewtonSolver {
    public:
        InnerNewtonSolver()
            : NewtonSolver()
        {
        }
        virtual ~InnerNewtonSolver() = default;

        virtual void c(const double value) { c_ = value; }
        virtual void e_b(const double value) { e_b_ = value; }
        virtual void t(const double value) { t_ = value; }
        virtual void m(const double value) { m_ = value; }

    protected:
        double line_search_lower_bound() const override
        {
            return std::min(
                NewtonSolver::line_search_lower_bound(), c_ * e_b_ / 10.0);
        }

        bool converged() override
        {
            switch (convergence_criteria) {
            case ConvergenceCriteria::ENERGY:
                spdlog::debug(
                    "solve={} iter={:d} step_energy={:g} tol={:g}", //
                    name(), iteration_number,
                    abs(gradient_free.dot(direction_free)),
                    std::max(e_b_, c_ * m_ / t_));
                return abs(direction_free.dot(gradient_free))
                    <= std::max(e_b_, c_ * m_ / t_);
            default:
                return NewtonSolver::converged();
            }
        }

        double c_;
        double e_b_;
        double t_;
        double m_;
    };

    std::unique_ptr<InnerNewtonSolver> inner_solver_ptr;
    BarrierProblem* problem_ptr;
    Eigen::VectorXd x0_i; ///< starting_point for each inner iteration
    int max_num_constraints;
    int num_outer_iterations;
};

} // namespace ipc::rigid
