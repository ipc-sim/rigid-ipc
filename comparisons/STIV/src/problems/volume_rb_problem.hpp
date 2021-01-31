#pragma once

#include <opt/optimization_problem.hpp>
#include <opt/volume_constraint.hpp>
#include <physics/rigid_body_problem.hpp>
#include <solvers/ncp_solver.hpp>

namespace ccd {

namespace opt {
    class VolumeRBProblem : public physics::RigidBodyProblem,
                            public virtual ConstrainedProblem {
    public:
        VolumeRBProblem(const std::string& name);

        void settings(const nlohmann::json& params) override;

        ////////////////////////////////////////////////////////////////
        /// Constrained Problem

        /// @brief eval_g evaluates constraints at point x
        virtual void compute_constraints(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& jac_gx,
            std::vector<Eigen::SparseMatrix<double>> hess_gx,
            bool compute_grad = true,
            bool compute_hess = true) override;

        using ConstrainedProblem::compute_constraints;

        virtual void compute_constraints_using_normal(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& jac_gx,
            std::vector<Eigen::SparseMatrix<double>> hess_gx,
            bool compute_grad = true,
            bool compute_hess = true) override;

        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }

        ////////////////////////////////////////////////////////////////
        /// ISimulationProblem
        opt::CollisionConstraint& constraint() override { return constraint_; }
        opt::OptimizationSolver& solver() override { return opt_solver_; }

    protected:
        void eval_jac_g_core(
            const Eigen::VectorXd& x,
            const EdgeEdgeImpacts&,
            Eigen::MatrixXd& jac_gx);

        opt::VolumeConstraint constraint_;
        opt::NCPSolver opt_solver_;
    };
} // namespace opt
} // namespace ccd
