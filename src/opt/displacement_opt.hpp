/**
 * Methods for optimizing the displacments with a non-linear interference volume
 * constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <opt/collision_constraint.hpp>
#include <opt/ncp_solver.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/solver.hpp>
#include <opt/volume_constraint.hpp>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * @brief Runs the optimization problem for the given method and initial
     * value
     * @param[in,out] problem : Optimization problem (possibly modified for
     *                          validation)
     * @param[in] U0          : Initial displacements
     * @param[in] settings    : Solver settings including tolerances and method
     * @return Optimization Result of the optimization
     * */
    OptimizationResults displacement_optimization(OptimizationProblem& problem,
        const Eigen::MatrixX2d& U0, SolverSettings& settings);

    class ParticlesDisplProblem : public OptimizationProblem {
    public:
        ParticlesDisplProblem();
        ~ParticlesDisplProblem() override;

        void initialize(const Eigen::MatrixX2d& V, const Eigen::MatrixX2i& E,
            const Eigen::MatrixX2d& U, CollisionConstraint& cstr);

        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;
        Eigen::MatrixX2d displacements;
        Eigen::MatrixXd u_;
        CollisionConstraint* constraint;

        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) override;
        Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        std::vector<Eigen::MatrixXd> eval_hessian_g(
            const Eigen::VectorXd& x) override;
        void initProblem();
    };

} // namespace opt
} // namespace ccd
