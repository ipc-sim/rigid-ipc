#pragma once

#include <opt/barrier_constraint.hpp>
#include <opt/particles_problem.hpp>

#include <physics/rigid_body_system.hpp>

namespace ccd {

namespace opt {

    class RigidBodyProblem2 : public OptimizationProblem {
    public:
        RigidBodyProblem2() {}
        virtual ~RigidBodyProblem2() override {}

        /// \brief initialize problem for new set of rigid bodies.
        void init(const std::vector<physics::RigidBody> rbs,
            opt::BarrierConstraint& constraint);

        /// \brief does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step);

        bool take_step(
            const Eigen::VectorXd& rb_positions, const double time_step);

        /// \brief update problem using current status of bodies.
        void update_constraint();

        double eval_f(const Eigen::VectorXd& sigma) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& sigma) override;
        Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& sigma) override;

        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& sigma) override
        {
            return Eigen::MatrixXd(eval_hessian_f_sparse(sigma));
        }

        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        bool has_barrier_constraint() override
        {
            return m_constraint_ptr->is_barrier();
        }
        double get_barrier_epsilon() override
        {
            return m_constraint_ptr->get_barrier_epsilon();
        }
        void set_barrier_epsilon(const double eps) override
        {
            return m_constraint_ptr->set_barrier_epsilon(eps);
        }

        physics::RigidBodyAssembler m_assembler;
        std::shared_ptr<opt::BarrierConstraint> m_constraint_ptr;

    protected:
        Eigen::MatrixXd m_x0; ///< vertices positions at begining of interval
        Eigen::VectorXd m_sigma1; ///< rigid body positions at end of interval
    };

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
            Eigen::VectorXd& g_uk,
            Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;
        ////////////////////////////////////////////////////////////////////////

        /// @brief Evaluate the intermediate callback.
        virtual bool eval_intermediate_callback(
            const Eigen::VectorXd& x) override;

    protected:
        virtual void init_num_vars() override;

        /// @brief pointer to the rigid bodies
        std::shared_ptr<ccd::physics::RigidBodySystem> model;
    };
} // namespace opt
} // namespace ccd
