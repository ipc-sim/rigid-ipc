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

    class ParticlesDisplProblem : public SimulationProblem {
    public:
        ParticlesDisplProblem();
        ParticlesDisplProblem(const std::string& name);
        ~ParticlesDisplProblem() override {}

        ////////////////////////////////////////////////////////////////////////
        /// SIMULATION

        void init(const nlohmann::json& params) override;
        nlohmann::json settings() const override;

        /// @brief  does a single simulation step. Returns true if there is a
        /// collision
        bool simulation_step(const double time_step) override;

        /// @brief moves status to given positions
        bool take_step(
            const Eigen::VectorXd& positions, const double time_step) override;

        /// \brief update optimization problem using current status.
        void update_constraint() override;

        /// --------------------------------------------------------------------------
        Eigen::MatrixXd vertices() const override { return vertices_; }
        /// @bried return candidate position of vertices, assuming no collisions
        Eigen::MatrixXd vertices_next(const double time_step) const override;
        Eigen::MatrixXd vertices_prev() const override
        {
            return vertices_prev_;
        }
        Eigen::MatrixXd vertices_collision() const override
        {
            return vertices_t1;
        }
        Eigen::MatrixXd dof_positions() const override { return vertices_; }

        /// \brief velocity of vertices at the END of the step
        /// as_delta = true, will return vertices / time_step;
        Eigen::MatrixXd velocities(
            const bool as_delta, const double time_step) const override;

        /// \brief collision forced applied to fix END position of vertices
        /// as_delta = true, will return Fc M^-1 / (time_step^2);
        Eigen::MatrixXd collision_force(
            const bool as_delta, const double time_step) const override;

        const Eigen::MatrixXi& edges() const override { return edges_; }

        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }
        const Eigen::MatrixXb& particle_dof_fixed() const override
        {
            return is_particle_dof_fixed;
        }
        const Eigen::VectorXd& gravity() const override { return gravity_; }

        bool detect_collisions(const Eigen::MatrixXd& q0,
            const Eigen::MatrixXd& q1,
            const CollisionCheck check_type) const;

        const opt::CollisionConstraint& constraint() override
        {
            return *constraint_ptr;
        }

        ////////////////////////////////////////////////////////////////////////
        /// CCD OPTIMIZATION PROBLEM
        ///
        /// @brief eval_f evaluates functional at point x
        virtual double eval_f(const Eigen::VectorXd& x) override;

        /// @brief eval_grad_f evaluates gradient of functional at point x
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        virtual Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        virtual void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override;
        ////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////
        // Constraint function and its derivatives.
        Eigen::MatrixXd update_g(const Eigen::VectorXd& x);

        /// @brief eval_g evaluates constraints at point x
        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        /// @brief eval_g evaluates constraints and jacobian at point x. Also
        /// returns list of active constraints (indices)
        virtual void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active) override;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        virtual Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        virtual void eval_jac_g(const Eigen::VectorXd& x,
            Eigen::SparseMatrix<double>& jac_gx) override;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        virtual std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        /// @brief eval_g_and_gdiff evaluates constraints, jacobian and hessian
        /// at point x
        virtual void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& g_uk,
            Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

        ///////////////////////////////////////////////////////////////////////
        /// FOR DEBUGGING
        ///
        /// creates sample points at xy coordinates
        void create_sample_points(const Eigen::MatrixXd& xy_points,
            Eigen::MatrixXd& sample_points) const override;

        void unflatten_dof(Eigen::MatrixXd& vec) const override
        {
            assert(vec.rows() % 2 == 0);
            vec.resize(vec.rows() / 2, 2);
        }

        ///////////////////////////////////////////////////////////////////////
        virtual bool has_barrier_constraint() override
        {
            return constraint_ptr->is_barrier();
        }
        virtual double get_barrier_epsilon() override
        {
            return constraint_ptr->get_barrier_epsilon();
        }
        virtual void set_barrier_epsilon(const double eps) override
        {
            return constraint_ptr->set_barrier_epsilon(eps);
        }
        ////////////////////////////////////////////////////////////////////////

        /// @brief Evaluate the intermediate callback.
        virtual bool eval_intermediate_callback(
            const Eigen::VectorXd& x) override;

        /// @brief Enable the line search mode. This functionality is not up to
        /// the child class.
        virtual void enable_line_search_mode(
            const Eigen::VectorXd& max_x) override;
        /// @brief Disable the line search mode.
        virtual void disable_line_search_mode() override;

        // ------------------------------------------------------------------------
        // World
        // ------------------------------------------------------------------------
        Eigen::MatrixXi edges_;
        Eigen::SparseMatrix<double> mass_matrix;
        Eigen::SparseMatrix<double> inv_mass_matrix;
        Eigen::MatrixXb is_particle_dof_fixed;
        Eigen::VectorXb is_dof_fixed_; ///> flattened version of above
        Eigen::VectorXd gravity_;

        // ------------------------------------------------------------------------
        // State
        // ------------------------------------------------------------------------
        Eigen::MatrixXd vertices_prev_; ///> vertices position at tau=0
        Eigen::MatrixXd vertices_;      ///> vertices position at tau=1
        Eigen::MatrixXd velocities_;    ///> velocities at tau=1
        Eigen::MatrixXd Fcollision;     ///> forced applied to solve collisions
        // ------------------------------------------------------------------------
        // Collision
        // ------------------------------------------------------------------------
        std::shared_ptr<opt::CollisionConstraint> constraint_ptr;
        intermediate_callback_func intermediate_callback;
        bool use_mass_matrix;
        double collision_eps;

        ////////////////////////////////////////////////////////////////////////

    protected:
        ///< 2D vertices positions at begining of interval
        Eigen::MatrixXd vertices_t0;
        ///< 2D vertices positions at end of interval
        Eigen::MatrixXd vertices_t1;
        ///< flatten vertices positions at begining of interval
        Eigen::MatrixXd vec_vertices_t0;
        ///< flatten vertices positions at end of interval
        Eigen::MatrixXd vec_vertices_t1;

        bool update_constraint_set;
        bool is_linesearch_active;
        ////////////////////////////////////////////////////////////////////////
    };

} // namespace physics
} // namespace ccd
