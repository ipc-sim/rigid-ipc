#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/volume_constraint.hpp>
#include <physics/particles_problem.hpp>
#include <solvers/ncp_solver.hpp>

namespace ccd {

namespace opt {
    class VolumeParticlesProblem : public physics::ParticlesProblem,
                                   public virtual opt::INCPProblem {
    public:
        VolumeParticlesProblem(const std::string& name);

        void settings(const nlohmann::json& params) override;

        ////////////////////////////////////////////////////////////////
        /// INCPProblem

        /// @brief eval_g evaluates constraints at point x
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        virtual void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active) override;

        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }

        ////////////////////////////////////////////////////////////////
        /// ISimulationProblem
        opt::CollisionConstraint& constraint() override { return constraint_; }
        opt::IStateOptimizationSolver& solver() override { return opt_solver_; }

    protected:
        opt::VolumeConstraint constraint_;
        opt::NCPSolver opt_solver_;
    };
} // namespace opt
} // namespace ccd
