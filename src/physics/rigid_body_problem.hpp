#pragma once

#include <memory> // shared_ptr

#include <physics/rigid_body_assembler.hpp>
#include <physics/simulation_problem.hpp>
#include <time_stepper/time_stepper.hpp>

namespace ccd {

namespace physics {

    class RigidBodyProblem : public virtual SimulationProblem {
    public:
        RigidBodyProblem();

        void init(const std::vector<RigidBody> rbs);

        virtual ~RigidBodyProblem() = default;

        static std::string problem_name() { return "rigid_body_problem"; }
        virtual std::string name() const override
        {
            return RigidBodyProblem::problem_name();
        }

        ////////////////////////////////////////////////////////////////////////
        /// Simulation Problem
        ////////////////////////////////////////////////////////////////////////

        virtual void settings(const nlohmann::json& params) override;
        nlohmann::json settings() const override;

        virtual nlohmann::json state() const override;
        void state(const nlohmann::json& s) override;

        /// @brief Compute the step but does not take it
        /// @returns true if there is a collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given configuration vector
        virtual bool
        take_step(const Eigen::VectorXd& x, const double time_step) override;

        /// @brief Check for intersections at the end of the time-step.
        virtual bool has_intersections() const override
        {
            return detect_intersections(this->poses_t1);
        }

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

        const Eigen::MatrixXi& edges() const override
        {
            return m_assembler.m_edges;
        }

        const Eigen::MatrixXi& faces() const override
        {
            return m_assembler.m_faces;
        }

        Eigen::MatrixXd velocities() const override
        {
            return m_assembler.world_velocities();
        }

        const Eigen::VectorXi& group_ids() const override
        {
            return m_assembler.group_ids();
        }

        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return m_assembler.is_dof_fixed;
        }

        Pose<double>
        rb_next_pose(const RigidBody& rb, const double time_step) const;

        /// Convert from dof expressed in distances to poses
        template <typename T>
        Eigen::VectorX<T> poses_to_dofs(const Poses<T>& poses) const
        {
            return Pose<T>::poses_to_dofs(poses);
        }

        /// Convert from poses to dof expressed in distances
        template <typename T>
        Poses<T> dofs_to_poses(const Eigen::VectorX<T>& dofs) const
        {
            return Pose<T>::dofs_to_poses(dofs, dim());
        }

        int dim() const override { return m_assembler.dim(); }

        virtual bool is_rb_problem() const override { return true; };

        // ------------------------------------------------------------------------
        // Settings
        // ------------------------------------------------------------------------
        double coefficient_restitution;
        Eigen::VectorXd gravity;
        double collision_eps;

        physics::RigidBodyAssembler m_assembler;

    protected:
        void solve_velocities();

        bool detect_collisions(
            const Poses<double>& poses_t0,
            const Poses<double>& poses_t1,
            const CollisionCheck check_type) const;

        bool detect_intersections(const Poses<double>& poses) const;

        void update_dof();

        /// @returns \f$x_0\f$: the starting point for the optimization.
        const Eigen::VectorXd& starting_point() const { return x0; }

        int num_vars_;
        /// Initial variable for optimization
        Eigen::VectorXd x0;

        /// Used during collision resolution
        /// Rigid body poses at start of time-step
        Poses<double> poses_t0;
        /// Rigid body poses at end of time-step
        Poses<double> poses_t1;

        /// Used for velocity restoration
        /// TODO: Replace this with the std::vector version
        ConcurrentImpacts original_impacts;

        std::shared_ptr<time_stepper::TimeStepper> m_time_stepper;

        /// Initial length of the bounding box diagonal
        double init_bbox_diagonal;
    };

    void assemble_hessian(
        const Eigen::SparseMatrix<double>& jac_xk_sigma,
        const std::vector<Eigen::SparseMatrix<double>>& hess_xk_sigma,
        const Eigen::MatrixXd& jac_g_uk,
        const std::vector<Eigen::SparseMatrix<double>>& hessian_g_uk,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian);

} // namespace physics
} // namespace ccd
