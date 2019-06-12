#pragma once

#include <opt/particles_problem.hpp>
#include <physics/rigid_body_system.hpp>

namespace ccd {

namespace opt {

    class RigidBodyProblem : public ParticlesDisplProblem {
    public:
        RigidBodyProblem();
        ~RigidBodyProblem() override;

        void initialize(
            ccd::physics::RigidBodySystem& rbs, CollisionConstraint& cstr);

        ////////////////////////////////////////////////////////////////////////
        // Objective function and its derivatives.
        /// @brief eval_f evaluates functional at point x
        virtual double eval_f(const Eigen::VectorXd& x) override;

        /// @brief eval_grad_f evaluates gradient of functional at point x
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        /// @brief eval_hessian_f evaluates hessian of functional at point x
        virtual Eigen::MatrixXd eval_hessian_f(
            const Eigen::VectorXd& x) override;

        /// @brief eval_hessian_f evaluates hessian of functional at point x
        virtual Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& x) override;
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Constraint function and its derivatives.
        /// @brief eval_g evaluates constraints at point x
        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        virtual Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        virtual std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        /// @brief eval_g_and_gdiff evaluates constraints, jacobian and hessian
        /// at point x
        virtual void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;
        ////////////////////////////////////////////////////////////////////////

    protected:
        virtual void init_num_vars() override;

        /// @brief pointer to the rigid bodies
        std::shared_ptr<ccd::physics::RigidBodySystem> model;
    };
} // namespace opt
} // namespace ccd
