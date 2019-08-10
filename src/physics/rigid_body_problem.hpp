#pragma once

#include <memory> // shared_ptr

#include <opt/collision_constraint.hpp>
#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>

namespace ccd {

namespace physics {

    class RigidBodyProblem : public SimulationProblem {
    public:
        RigidBodyProblem(const std::string& name);
        RigidBodyProblem();

        virtual ~RigidBodyProblem() override {}

        ////////////////////////////////////////////////////////////////////////
        // SIMULATION
        ////////////////////////////////////////////////////////////////////////

        /// \brief initialize problem for new set of rigid bodies.
        void init(const nlohmann::json& params) override;
        void init(
            const std::vector<RigidBody> rbs, const std::string& constraint);
        nlohmann::json settings() const override;
        nlohmann::json state() const override;
        void state(const nlohmann::json& s) override;

        /// Update state of system
        ///-----------------------------------------------------------------

        /// \brief does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given positions
        bool take_step(const Eigen::VectorXd& rb_positions,
            const double time_step) override;

        /// \brief update problem using current status of bodies.
        void update_constraint() override;

        /// GET state of system
        ///-----------------------------------------------------------------
        bool detect_collisions(const Eigen::MatrixXd& q0,
            const Eigen::MatrixXd& q1,
            const CollisionCheck check_type);

        /// \brief returns world vertices at the END of NEXT step.
        /// Returns the positions of the next step assuming no collisions forces
        Eigen::MatrixXd vertices_next(const double time_step) const override;

        /// \brief returns world vertices at the END of step (current)
        Eigen::MatrixXd vertices() const override
        {
            return m_assembler.world_vertices_t1();
        }
        /// \brief returns world vertices at the BEGINNING of step (current)
        Eigen::MatrixXd vertices_prev() const override
        {
            return m_assembler.world_vertices_t0();
        }
        /// \brief returns world vertices that we are solving collisions for
        Eigen::MatrixXd vertices_collision() const override
        {
            return vertices_q1;
        }

        Eigen::MatrixXd dof_positions() const override
        {
            Eigen::MatrixXd p = m_assembler.rb_positions_t1();
            unflatten_dof(p);
            return p;
        }

        Eigen::Vector3d rb_position_next(
            const RigidBody& rb, const double time_step) const;

        /// \brief velocity of vertices at the END of the step
        /// as_delta = true, will return vertices / time_step;
        Eigen::MatrixXd velocities(
            const bool as_delta, const double time_step) const override;

        /// \brief collision forced applied to fix END position of vertices
        /// as_delta = true, will return Fc M^-1 / (time_step^2);
        Eigen::MatrixXd collision_force(
            const bool as_delta, const double time_step) const override;

        void unflatten_dof(Eigen::MatrixXd& vec) const override
        {
            // order is x,y,z,x,y,z,x,y,z,
            assert(vec.rows() % 3 == 0);
            vec.resize(3, vec.rows() / 3);
            vec.transposeInPlace();
        }

        const Eigen::MatrixXi& edges() const override
        {
            return m_assembler.m_edges;
        }

        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return m_assembler.is_dof_fixed;
        }
        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }

        const Eigen::VectorXd& gravity() const override { return gravity_; }

        const opt::CollisionConstraint& constraint() override
        {
            return *m_constraint_ptr;
        }

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

        void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active) override;

        void eval_jac_g(const Eigen::VectorXd& x,
            Eigen::SparseMatrix<double>& jac_gx) override;

        Eigen::VectorXd eval_g(const Eigen::VectorXd& x,
            const bool update_constraint_set) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x,
            const bool update_constraint_set) override;
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x,
            const bool update_constraint_set) override;

        Eigen::MatrixXd update_g(const Eigen::VectorXd& gamma);
        Eigen::MatrixXd update_g(
            const Eigen::VectorXd& gamma, const bool update_constraint_set);

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
        bool m_update_constraint_set;
        double coefficient_restitution;
        Eigen::VectorXd gravity_;
        double collision_eps;

    protected:
        void solve_velocities();
        /// Used during collision resolution
        ///< vertices positions at begining of interval
        Eigen::MatrixXd vertices_t0;
        ///< vertices positions at end of interval
        Eigen::MatrixXd vertices_q1;
        ///< rigid body positions at end of interval
        Eigen::VectorXd rb_positions_t1;

        /// Used for velocity restoration
        EdgeVertexImpacts original_ev_impacts;

        /// Used for visualization and debugging
        Eigen::MatrixXd Fcollision; ///< forces used to resolve collisions
    };

} // namespace physics
} // namespace ccd
