#pragma once

#include <memory> // shared_ptr

#include <opt/collision_constraint.hpp>
#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>

namespace ccd {

namespace physics {

    class RigidBodyProblem : public SimulationProblem {
    public:
        RigidBodyProblem();

        virtual ~RigidBodyProblem() override {}

        ////////////////////////////////////////////////////////////////////////
        // SIMULATION
        ////////////////////////////////////////////////////////////////////////
        bool validate_params(const nlohmann::json& params);

        /// \brief initialize problem for new set of rigid bodies.
        void init(const nlohmann::json& params) override;
        void init(
            const std::vector<RigidBody> rbs, const std::string& constraint);

        /// \brief does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        bool take_step(const Eigen::VectorXd& rb_positions,
            const double time_step) override;

        bool detect_collisions(const Eigen::MatrixXd& q0,
            const Eigen::MatrixXd& q1,
            const CollisionCheck check_type);

        /// \brief returns world vertices at the END of NEXT step.
        /// Returns the positions of the next step assuming no collisions forces
        Eigen::MatrixXd vertices_next(const double time_step) override;
        Eigen::Vector3d rb_position_next(
            const RigidBody& rb, const double time_step) const;

        /// \brief update problem using current status of bodies.
        void update_constraint() override;

        Eigen::MatrixXd velocities(
            const bool as_delta, const double time_step) override;

        Eigen::MatrixXd collision_force(
            const bool as_delta, const double time_step) override;

        const opt::CollisionConstraint& constraint() override
        {
            return *m_constraint_ptr;
        }

        /// \brief returns world vertices at the END of step (current)
        Eigen::MatrixXd vertices() override
        {
            return m_assembler.world_vertices_t1();
        }
        /// \brief returns world vertices at the BEGINNING of step (current)
        Eigen::MatrixXd vertices_prev() override
        {
            return m_assembler.world_vertices_t0();
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
        const Eigen::VectorXd& gravity() override { return gravity_; }

        ////////////////////////////////////////////////////////////////////////
        // Objective function and its derivatives.
        ////////////////////////////////////////////////////////////////////////
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
        ////////////////////////////////////////////////////////////////////////
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
        /// FOR DEBUGGING
        ///
        /// creates sample points at xy coordinates
        void create_sample_points(const Eigen::MatrixXd& xy_points,
            Eigen::MatrixXd& sample_points) override;

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
        Eigen::VectorXd gravity_;
        double collision_eps;

    protected:
        Eigen::MatrixXd m_q0; ///< vertices positions at begining of interval
        Eigen::MatrixXd m_q1; ///< vertices positions at end of interval
        Eigen::MatrixXd m_Fcollision; ///< forces used to resolve collisions
        Eigen::VectorXd m_sigma1; ///< rigid body positions at end of interval
    };

} // namespace physics
} // namespace ccd
