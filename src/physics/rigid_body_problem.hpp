#pragma once

#include <memory> // shared_ptr

#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>

namespace ccd {

namespace physics {

    class RigidBodyProblem : public virtual ISimulationProblem,
                             public virtual opt::IUnconstraintedProblem {
    public:
        RigidBodyProblem(const std::string& name);
        RigidBodyProblem();

        virtual ~RigidBodyProblem() override = default;

        ////////////////////////////////////////////////////////////////////////
        /// I-SIMULATION
        ////////////////////////////////////////////////////////////////////////

        // virtual functions to implement by child classes
        // ----------------------------------------------------
        // virtual opt::CollisionConstraint& constraint() = 0;
        // virtual opt::IStateOptimizationSolver& solver() = 0;

        std::string name() override { return name_; }
        virtual void settings(const nlohmann::json& params) override;
        nlohmann::json settings() const override;
        virtual nlohmann::json state() const override;
        void state(const nlohmann::json& s) override;

        /// \brief does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given configuration vector
        virtual bool take_step(const Eigen::VectorXd& sigma,
            const double time_step) override;

        /// \brief update problem using current status of bodies.
        void update_constraint() override;
        opt::OptimizationResults solve_constraints() override;
        void init_solve() override;
        opt::OptimizationResults step_solve() override;

        /// \brief returns world vertices at the END of step (current)
        Eigen::MatrixXd vertices() const override
        {
            return m_assembler.world_vertices_t1();
        }

        Eigen::MatrixXd velocities() const override {
            return  m_assembler.world_velocities();
        }
        const Eigen::MatrixXi& edges() const override
        {
            return m_assembler.m_edges;
        }

        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return m_assembler.is_dof_fixed;
        }
        void init(const std::vector<RigidBody> rbs);

        Eigen::Vector3d rb_position_next(
            const RigidBody& rb, const double time_step) const;

        ////////////////////////////////////////////////////////////////////////
        /// IUnconstraintedProblem
        ////////////////////////////////////////////////////////////////////////

        /// @brief eval_f evaluates functional at point x
        double eval_f(const Eigen::VectorXd& sigma) override;

        /// @brief eval_grad_f evaluates gradient of functional at point x
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& sigma) override;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& sigma) override;

        const int& num_vars() override { return num_vars_; }
        const Eigen::VectorXd& starting_point() override { return x0; }

        // ------------------------------------------------------------------------
        // Settings
        // ------------------------------------------------------------------------
        double coefficient_restitution;
        Eigen::VectorXd gravity_;
        double collision_eps;

        physics::RigidBodyAssembler m_assembler;

    protected:
        void solve_velocities();

        bool detect_collisions(const Eigen::MatrixXd& q0,
            const Eigen::MatrixXd& q1,
            const CollisionCheck check_type);

        int num_vars_;
        Eigen::VectorXd x0;

        /// Used during collision resolution
        ///< vertices positions at begining of interval
        Eigen::MatrixXd vertices_t0;
        ///< vertices positions at end of interval
        Eigen::MatrixXd vertices_q1;
        ///< rigid body positions at end of interval
        Eigen::VectorXd sigma_t1;

        /// Used for velocity restoration
        EdgeVertexImpacts original_ev_impacts;

        /// Used for visualization and debugging
        Eigen::MatrixXd Fcollision; ///< forces used to resolve collisions

        std::string name_;
    };

    void assemble_hessian(const Eigen::SparseMatrix<double>& jac_xk_sigma,
        const std::vector<Eigen::SparseMatrix<double>>& hess_xk_sigma,
        const Eigen::MatrixXd& jac_g_uk,
        const std::vector<Eigen::SparseMatrix<double>>& hessian_g_uk,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian);

} // namespace physics
} // namespace ccd
