/**
 * Methods for optimizing the displacements with a non-linear interference
 * volume constraint.
 */

#pragma once

#include <Eigen/Core>
#include <memory> // shared_ptr

#include <ccd/collision_detection.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/volume_constraint.hpp>
#include <solvers/optimization_solver.hpp>

#include <physics/simulation_problem.hpp>

#include <utils/not_implemented_error.hpp>

namespace ccd {

namespace physics {

    typedef std::function<bool(const Eigen::VectorXd&, const Eigen::MatrixX2d&)>
        intermediate_callback_func;

    class ParticlesProblem : public virtual ISimulationProblem,
                             public virtual opt::IUnconstraintedProblem {
    public:
        ParticlesProblem();
        ParticlesProblem(const std::string& name);
        ~ParticlesProblem() override {}

        ////////////////////////////////////////////////////////////////////////
        /// I-SIMULATION
        ////////////////////////////////////////////////////////////////////////

        // virtual functions to implement by child classes
        // ----------------------------------------------------
        // virtual opt::CollisionConstraint& constraint() = 0;
        // virtual opt::IStateOptimizationSolver& solver() = 0;

        std::string name() override { return name_; }
        void settings(const nlohmann::json& params) override;
        nlohmann::json settings() const override;
        nlohmann::json state() const override;
        void state(const nlohmann::json& s) override;

        /// @brief  does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given positions
        bool take_step(
            const Eigen::VectorXd& positions, const double time_step) override;

        /// \brief update problem using current status of bodies.
        void update_constraint() override;
        opt::OptimizationResults solve_constraints() override;
        void init_solve() override;
        opt::OptimizationResults step_solve() override;

        Eigen::MatrixXd vertices() const override { return vertices_; }
        Eigen::MatrixXd velocities() const override { return velocities_; }

        const Eigen::MatrixXi& edges() const override { return edges_; }

        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return is_particle_dof_fixed;
        }

        Eigen::VectorXi group_id() const override { return group_ids_; }

        ////////////////////////////////////////////////////////////////////////
        /// IUnconstraintedProblem

        /// @brief eval_f evaluates functional at point x
        virtual double eval_f(const Eigen::VectorXd& x) override;

        /// @brief eval_grad_f evaluates gradient of functional at point x
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        virtual Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        const int& num_vars() override { return num_vars_; }
        const Eigen::VectorXd& starting_point() override { return x0; }

        // ------------------------------------------------------------------------
        // World
        // ------------------------------------------------------------------------
        Eigen::MatrixXi edges_;
        Eigen::VectorXi group_ids_;
        Eigen::SparseMatrix<double> mass_matrix;
        Eigen::MatrixXb is_particle_dof_fixed;
        Eigen::VectorXb is_dof_fixed_; ///> flattened version of above
        Eigen::VectorXd gravity_;

        // ------------------------------------------------------------------------
        // State
        // ------------------------------------------------------------------------
        Eigen::MatrixXd vertices_prev_; ///> vertices position at tau=0
        Eigen::MatrixXd vertices_;      ///> vertices position at tau=1
        Eigen::MatrixXd velocities_;    ///> velocities at tau=1

        // ------------------------------------------------------------------------
        // Collision
        // ------------------------------------------------------------------------

        double collision_eps;

        ////////////////////////////////////////////////////////////////////////

    protected:
        Eigen::MatrixXd vertices_next(const double time_step) const;

        bool detect_collisions(const Eigen::MatrixXd& q0,
            const Eigen::MatrixXd& q1,
            const CollisionCheck check_type);

        int num_vars_;
        int num_constraints_;
        Eigen::VectorXd x0;

        ///< 2D vertices positions at begining of interval
        Eigen::MatrixXd vertices_t0;
        ///< 2D vertices positions at end of interval
        Eigen::MatrixXd vertices_t1;
        ///< flatten vertices positions at begining of interval
        Eigen::MatrixXd vec_vertices_t0;
        ///< flatten vertices positions at end of interval
        Eigen::MatrixXd vec_vertices_t1;

        std::string name_;
        ////////////////////////////////////////////////////////////////////////
    };

} // namespace physics
} // namespace ccd
