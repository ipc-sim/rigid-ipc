#pragma once

#include <autodiff/autodiff_types.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <physics/rigid_body_problem.hpp>
#include <solvers/barrier_solver.hpp>

namespace ccd {

namespace opt {

    struct RB2Candidate {
        int vertex_body_id;
        int edge_body_id;
        int vertex_local_id;
        int edge0_local_id;
        int edge1_local_id;
    };

    class DistanceBarrierRBProblem
        : public physics::RigidBodyProblem,
          public virtual opt::IBarrierGeneralProblem {
    public:
        DistanceBarrierRBProblem(const std::string& name);

        void settings(const nlohmann::json& params) override;

        ////////////////////////////////////////////////////////////////
        /// IConstrainedProblem

        /// @brief eval_g evaluates constraints at point x
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        const int& num_constraints() override
        {
            return constraint_.m_num_constraints;
        }

        ////////////////////////////////////////////////////////////
        /// IBarrierProblem

        Eigen::VectorXd eval_g_set(
            const Eigen::VectorXd& x, const CstrSetFlag flag) override;

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
            return constraint_.m_barrier_epsilon;
        }

        void set_barrier_epsilon(const double eps) override
        {
            constraint_.m_barrier_epsilon = eps;
        }

        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }

        opt::CollisionConstraint& constraint() override { return constraint_; }
        opt::IStateOptimizationSolver& solver() override { return opt_solver_; }

    protected:
        void update_constraints(
            const Eigen::MatrixXd& xk, const CstrSetFlag flag);

        void extract_local_system(
            const EdgeVertexCandidate& c, RB2Candidate& rbc);

        template <typename T>
        T distance_barrier(
            const Eigen::VectorXd& sigma, const RB2Candidate& rbc);

        Eigen::MatrixXd eval_jac_g_core(const Eigen::VectorXd& sigma);
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g_core(
            const Eigen::VectorXd& sigma);

        opt::DistanceBarrierConstraint constraint_;
        opt::BarrierSolver opt_solver_;
    };

} // namespace opt
} // namespace ccd
