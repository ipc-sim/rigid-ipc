#pragma once

#include <opt/collision_constraint.hpp>
#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>

namespace ccd {

namespace physics {

    class RigidBodyProblem : public SimulationProblem {
    public:
        RigidBodyProblem();

        virtual ~RigidBodyProblem() override {}

        /// \brief initialize problem for new set of rigid bodies.
        void init(const nlohmann::json& params) override;

        void init(
            const std::vector<RigidBody> rbs, const std::string& constraint);

        /// \brief does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        bool take_step(const Eigen::VectorXd& rb_positions,
            const double time_step) override;

        /// \brief update problem using current status of bodies.
        void update_constraint() override;

        const opt::CollisionConstraint& constraint() override
        {
            return *m_constraint_ptr;
        }

        Eigen::MatrixXd vertices() override
        {
            return m_assembler.world_vertices_t1();
        }
        const Eigen::MatrixXi& edges() override { return m_assembler.m_edges; }
        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }
        const Eigen::MatrixXb& particle_dof_fixed() override
        {
            return m_assembler.is_dof_fixed;
        }

        ////////////////////////////////////////////////////////////////////////
        // Objective function and its derivatives.
        double eval_f(const Eigen::VectorXd& sigma) override;

        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& sigma) override;

        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& sigma) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;

        /// \brief functional problem using chain rule
        double eval_f_chain(const Eigen::VectorXd& sigma);

        Eigen::VectorXd eval_grad_f_chain(const Eigen::VectorXd& sigma);

        Eigen::SparseMatrix<double> eval_hessian_f_chain(
            const Eigen::VectorXd& sigma);

        ////////////////////////////////////////////////////////////////////////
        // Constraint function and its derivatives.
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& gx_hessian) override;

        Eigen::MatrixXd update_g(const Eigen::VectorXd& gamma);

        ///////////////////////////////////////////////////////////////////////
        /// BARRIER SPECIFIC
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
        std::shared_ptr<opt::CollisionConstraint> m_constraint_ptr;
        bool use_chain_functional;
        bool update_constraint_set;

    protected:
        Eigen::MatrixXd m_q0; ///< vertices positions at begining of interval
        Eigen::MatrixXd m_q1; ///< vertices positions at end of interval
        Eigen::VectorXd m_sigma1; ///< rigid body positions at end of interval
    };

} // namespace physics
} // namespace ccd
