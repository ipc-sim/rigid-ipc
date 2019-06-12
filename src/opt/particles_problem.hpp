/**
 * Methods for optimizing the displacements with a non-linear interference
 * volume constraint.
 */

#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

#include <opt/collision_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <opt/volume_constraint.hpp>

#include <solvers/optimization_solver.hpp>

namespace ccd {

namespace opt {

    typedef std::function<bool(const Eigen::VectorXd&, const Eigen::MatrixX2d&)>
        intermediate_callback_func;

    class ParticlesDisplProblem : public OptimizationProblem {
    public:
        ParticlesDisplProblem();
        virtual ~ParticlesDisplProblem() override;

        void initialize(const Eigen::MatrixX2d& V, const Eigen::MatrixX2i& E,
            const Eigen::MatrixX2d& U, CollisionConstraint& cstr);

        /// \brief eval_f evaluates functional at point x
        double eval_f(const Eigen::VectorXd& x) override;

        /// \brief eval_grad_f evaluates gradient of functional at point x
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;

        /// \brief eval_hessian_f evaluates hessian of functional at point x
        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) override;
        Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& x) override;

        /// \brief eval_g evaluates constraints at point x
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

        /// \brief eval_g evaluates constraints and jacobian at point x. Also
        /// returns list of active constraints (indices)
        void eval_g(const Eigen::VectorXd& x, Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active) override;

        /// \brief eval_g_and_gdiff evaluates constraints, jacobian and hessian
        /// at point x
        void eval_g_and_gdiff(const Eigen::VectorXd& x, Eigen::VectorXd& g_uk,
            Eigen::MatrixXd& g_uk_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian) override;

        void enable_line_search_mode(const Eigen::VectorXd& max_x) override;
        void disable_line_search_mode() override;

        /// \brief eval_jac_g evaluates constraints jacobian at point x
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        void eval_jac_g(const Eigen::VectorXd& x,
            Eigen::SparseMatrix<double>& jac_gx) override;

        // \brief eval_hessian_g evaluates constraints hessian at point x
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        bool eval_intermediate_callback(const Eigen::VectorXd& x) override;

        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;
        Eigen::MatrixX2d displacements;
        Eigen::MatrixXd u_;
        CollisionConstraint* constraint;
        intermediate_callback_func intermediate_callback;
        bool use_mass_matrix;

    private:
        bool is_collision_set_frozen;

        void initProblem();

        Eigen::SparseMatrix<double> mass_matrix;
        void init_mass_matrix();
    };

} // namespace opt
} // namespace ccd
