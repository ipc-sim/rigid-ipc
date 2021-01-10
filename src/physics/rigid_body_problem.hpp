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

        virtual void init(const std::vector<RigidBody>& rbs);

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

        virtual double timestep() const override { return m_timestep; }
        virtual void timestep(double timestep) override
        {
            m_timestep = timestep;
        }

        /// Update problem using current status of bodies.
        virtual void update_constraints();
        /// Initialize the solver.
        virtual void init_solve();
        /// Use the solver to solve this problem.
        virtual opt::OptimizationResults solve_constraints();
        /// Take a single solver step.
        virtual opt::OptimizationResults step_solve();

        /// World vertices at the END of step (current).
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

        const Eigen::MatrixXb& vertex_dof_fixed() const override
        {
            return m_assembler.is_dof_fixed;
        }

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
        long num_vertices() const override
        {
            return m_assembler.num_vertices();
        }
        long num_edges() const override { return m_assembler.num_edges(); }
        long num_faces() const override { return m_assembler.num_faces(); }
        size_t num_bodies() const override { return m_assembler.num_bodies(); }

        virtual bool is_rb_problem() const override { return true; };

        // --------------------------------------------------------------------
        // Settings
        // --------------------------------------------------------------------
        double coefficient_restitution; ///< Coefficent of resitution
        double coefficient_friction;    ///< Coefficent of friction
        Eigen::VectorX3d gravity;       ///< Acceleration due to gravity
        double collision_eps; ///< Scale trajectory for early collision

        physics::RigidBodyAssembler m_assembler;

    protected:
        /// Moves status to given configuration vector.
        virtual bool take_step(const Eigen::VectorXd& x);

        /// Detect collisions between poses_t0 and poses_t1.
        bool detect_collisions(
            const Poses<double>& poses_t0,
            const Poses<double>& poses_t1,
            const CollisionCheck check_type) const;

        /// Detect intersections between rigid bodies with given poses.
        bool detect_intersections(const Poses<double>& poses) const;

        virtual void update_dof();

        /// @returns \f$x_0\f$: the starting point for the optimization.
        const Eigen::VectorXd& starting_point() const { return x0; }

        double m_timestep; ///< The time-step size

        Eigen::VectorXd x0; ///< Initial variable for optimization
        int num_vars_;      ///< The number of variables

        // Used during collision resolution
        Poses<double> poses_t0; ///< Rigid body poses at start of time-step

        /// Initial length of the bounding box diagonal
        double init_bbox_diagonal;
    };

} // namespace physics
} // namespace ccd
