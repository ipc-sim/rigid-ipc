#pragma once

#include <memory> // shared_ptr

#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>
#include <time_stepper/time_stepper.hpp>

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

        std::string name() override { return name_; }
        virtual void settings(const nlohmann::json& params) override;
        nlohmann::json settings() const override;
        virtual nlohmann::json state() const override;
        void state(const nlohmann::json& s) override;

        /// @brief take a single simulation step
        /// @returns true if there is a collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given configuration vector
        virtual bool
        take_step(const Eigen::VectorXd& dof, const double time_step) override;

        /// @brief update problem using current status of bodies.
        void update_constraint() override;
        opt::OptimizationResults solve_constraints() override;
        void init_solve() override;
        opt::OptimizationResults step_solve() override;

        /// @brief returns world vertices at the END of step (current)
        Eigen::MatrixXd vertices() const override
        {
            return m_assembler.world_vertices_t1();
        }

        Eigen::MatrixXd velocities() const override
        {
            return m_assembler.world_velocities();
        }

        const Eigen::MatrixXi& edges() const override
        {
            return m_assembler.m_edges;
        }

        const Eigen::MatrixXi& faces() const override
        {
            return m_assembler.m_faces;
        }

        Eigen::VectorXi group_ids() const override
        {
            return m_assembler.group_ids();
        }

        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return m_assembler.is_dof_fixed;
        }
        void init(const std::vector<RigidBody> rbs);

        Pose<double>
        rb_next_pose(const RigidBody& rb, const double time_step) const;

        /// Convert from dof expressed in distances to poses
        template <typename T>
        Eigen::VectorX<T> poses_to_dofs(const Poses<T>& poses) const
        {
            return m_assembler.m_pose_to_dof.cast<T>()
                * Pose<T>::poses_to_dofs(poses);
        }

        /// Convert from poses to dof expressed in distances
        template <typename T>
        Poses<T> dofs_to_poses(const Eigen::VectorX<T>& dofs) const
        {
            return Pose<T>::dofs_to_poses(
                m_assembler.m_dof_to_pose.cast<T>() * dofs, dim());
        }

        ////////////////////////////////////////////////////////////////////////
        /// IUnconstraintedProblem
        ////////////////////////////////////////////////////////////////////////

        /// @brief eval_f evaluates functional at point x
        double eval_f(const Eigen::VectorXd& sigma) override;

        /// @brief eval_grad_f evaluates gradient of functional at point x
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& sigma) override;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        Eigen::SparseMatrix<double>
        eval_hessian_f(const Eigen::VectorXd& sigma) override;

        const int& num_vars() override { return num_vars_; }
        const Eigen::VectorXd& starting_point() override { return x0; }
        int dim() const { return m_assembler.dim(); }

        // ------------------------------------------------------------------------
        // Settings
        // ------------------------------------------------------------------------
        double coefficient_restitution;
        Eigen::VectorXd gravity_;
        double collision_eps;

        physics::RigidBodyAssembler m_assembler;

    protected:
        void solve_velocities();

        bool detect_collisions(
            const Poses<double>& poses_t0,
            const Poses<double>& poses_t1,
            const CollisionCheck check_type) const;

        void update_dof();

        int num_vars_;
        Eigen::VectorXd x0;

        /// Used during collision resolution
        /// Rigid body poses at start of time-step
        Poses<double> poses_t0;
        /// Rigid body poses at end of time-step
        Poses<double> poses_t1;

        /// Used for velocity restoration
        /// TODO: Replace this with the std::vector version
        ConcurrentImpacts original_impacts;

        std::string name_;

        std::shared_ptr<time_stepper::TimeStepper> m_time_stepper;
    };

    void assemble_hessian(
        const Eigen::SparseMatrix<double>& jac_xk_sigma,
        const std::vector<Eigen::SparseMatrix<double>>& hess_xk_sigma,
        const Eigen::MatrixXd& jac_g_uk,
        const std::vector<Eigen::SparseMatrix<double>>& hessian_g_uk,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian);

} // namespace physics
} // namespace ccd
