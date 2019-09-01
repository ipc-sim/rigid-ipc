#pragma once

#include <autodiff/autodiff_types.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <physics/particles_problem.hpp>
#include <solvers/barrier_solver.hpp>

namespace ccd {
namespace opt {
    class DistanceBarrierParticleProblem
        : public physics::ParticlesProblem,
          public virtual opt::IBarrierGeneralProblem {
    public:
        DistanceBarrierParticleProblem(const std::string& name);

        void settings(const nlohmann::json& params) override;

        ////////////////////////////////////////////////////////////////
        /// IConstrainedProblem
        ////////////////////////////////////////////////////////////////

        /// @brief eval_g evaluates constraints at point x
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        const int& num_constraints() override
        {
            throw NotImplementedError(
                "DistanceBarrierParticleProblem non-const num_constraints");
        }
        ////////////////////////////////////////////////////////////
        /// IBarrierProblem
        ////////////////////////////////////////////////////////////////
        bool has_collisions(const Eigen::VectorXd& sigma_i,
            const Eigen::VectorXd& sigma_j) const override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian) override;

        void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& gx_hessian) override;

        double get_barrier_epsilon() override
        {
            return constraint_.get_barrier_epsilon();
        }

        void set_barrier_epsilon(const double eps) override
        {
            constraint_.set_barrier_epsilon(eps);
        }

        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }

        opt::CollisionConstraint& constraint() override { return constraint_; }
        opt::IStateOptimizationSolver& solver() override { return opt_solver_; }

    protected:
        opt::DistanceBarrierConstraint constraint_;
        opt::BarrierSolver opt_solver_;
    };
} // namespace opt
} // namespace ccd
