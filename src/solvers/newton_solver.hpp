#pragma once

#include <Eigen/Core>

#include <constants.hpp>
#include <solvers/optimization_solver.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

enum ConvergenceCriteria {
    VELOCITY, ///< Change in position of a newton direction
    ENERGY    ///< Change in energy of a newton direction
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    ConvergenceCriteria, { { VELOCITY, "velocity" }, { ENERGY, "energy" } });

class NewtonSolver : public virtual OptimizationSolver {
public:
    NewtonSolver();
    virtual ~NewtonSolver() = default;

    /// Initialize the state of the solver using the settings saved in JSON
    virtual void settings(const nlohmann::json& params) override;
    /// Export the state of the solver using the settings saved in JSON
    virtual nlohmann::json settings() const override;

    /// An identifier for the solver class
    static std::string solver_name() { return "newton_solver"; }
    /// An identifier for this solver
    virtual std::string name() const override
    {
        return NewtonSolver::solver_name();
    }

    /// Initialize the solver with a problem to solve
    virtual void set_problem(OptimizationProblem& problem) override
    {
        this->problem_ptr = &problem;
    }

    /// Initialize the solver state for a new solve
    virtual void init_solve(const Eigen::VectorXd& x0) override;

    /**
     * @brief Perform Newton's Method to minimize the objective, \f$f(x)\f$,
     * of the problem unconstrained.
     *
     * @param[in] problem  The optimization problem to minimize
     *                     unconstrained.
     *
     * @return The results of the optimization including the minimizer,
     * minimum, and if the optimization was successful.
     */
    virtual OptimizationResults solve(const Eigen::VectorXd& x0) override;

    /// Perform a single step of solving the optimization problem
    virtual OptimizationResults step_solve() override
    {
        throw NotImplementedError(
            "Taking a single newton step is not implemented!");
    };

    /**
     * @brief Solve for the Newton direction
     *        (\f$\Delta x = -H^{-1} \nabla f \f$).
     *
     * @param[in]  gradient  Gradient of the objective function.
     * @param[in]  hessian   Hessian of the objective function.
     * @param[out] delta_x   Output newton direction.
     * @param[in]  make_psd  If delta_x is not a descent direction, then
     * make the hessian positive semi-definite.
     *
     * @return Returns true if the solve was successful.
     */
    virtual bool compute_direction(
        const Eigen::VectorXd& gradient,
        const Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& delta_x,
        bool make_psd = false);

    virtual bool compute_regularized_direction(
        double& fx,
        Eigen::VectorXd& gradient,
        Eigen::SparseMatrix<double>& hessian,
        Eigen::VectorXd& delta_x,
        double& coeff);

    virtual std::string stats_string() const override;
    virtual nlohmann::json stats() const override;

    int max_iterations;

protected:
    virtual bool converged();

    virtual void post_step_update();

    virtual bool line_search(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const double fx,
        const Eigen::VectorXd& grad_fx,
        double& step_length);

    virtual double line_search_lower_bound() const
    {
        return m_line_search_lower_bound;
    }

    /// @brief Pointer to the problem to solve.
    OptimizationProblem* problem_ptr;

    int iteration_number; ///< @brief The current iteration number.
    ConvergenceCriteria convergence_criteria;

    double m_line_search_lower_bound; ///< @brief Line search lower bound

    double energy_conv_tol;        ///< @brief Energy convergence tolerance
    double velocity_conv_tol;      ///< @brief Velocity convergence tolerance
    bool is_velocity_conv_tol_abs; ///< @brief Absolute velocity tol
    bool is_energy_converged;

    // State variables
    Eigen::VectorXd x, x_prev;
    Eigen::VectorXd gradient, gradient_free;
    Eigen::VectorXd direction, direction_free;
    Eigen::VectorXd grad_direction; ///< Gradient with fixed DoF set to zero
    Eigen::SparseMatrix<double> hessian, hessian_free;

private:
    void reset_stats();

    size_t num_fx = 0;
    size_t num_grad_fx = 0;
    size_t num_hessian_fx = 0;
    size_t num_collision_check = 0;
    size_t ls_iterations = 0;
    size_t newton_iterations = 0;
    size_t num_newton_ls_fails = 0;
    size_t num_grad_ls_fails = 0;
    size_t regularization_iterations = 0;
};

/**
 * @brief Make the matrix positive definite (\f$x^T A x > 0\$).
 *
 * @param A The matrix to make positive definite.
 *
 * @return The scale of the update to the diagonal.
 */
double make_matrix_positive_definite(Eigen::SparseMatrix<double>& A);

/**
 * @brief Log values along a search direction.
 *
 * @param[in] x                 Starting point for the line search.
 * @param[in] dir               Direction to search along.
 * @param[in] f_and_gradf       Function of x to sample with gradient.
 */
void sample_search_direction(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dir,
    const std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)>&
        f_and_gradf,
    double max_step);

} // namespace ipc::rigid
