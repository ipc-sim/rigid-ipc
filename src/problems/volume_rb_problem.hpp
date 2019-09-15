#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/volume_constraint.hpp>
#include <physics/rigid_body_problem.hpp>
#include <solvers/ncp_solver.hpp>

namespace ccd {

namespace opt {
    class VolumeRBProblem : public physics::RigidBodyProblem,
                            public virtual opt::INCPProblem {
    public:
        VolumeRBProblem(const std::string& name);

        void settings(const nlohmann::json& params) override;

        ////////////////////////////////////////////////////////////////
        /// INCPProblem

        /// @brief eval_g evaluates constraints at point x
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        virtual void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian) override;

        virtual void eval_g_normal(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian) override;


        void eval_jac_g_core(const Eigen::VectorXd& x,
            const EdgeEdgeImpacts&,
            Eigen::MatrixXd& jac_gx);

        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }

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
