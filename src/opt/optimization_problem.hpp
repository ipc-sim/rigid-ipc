#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    /// @brief default value for no upper bound
    static const double NO_UPPER_BOUND = 2e19;
    /// @brief default value for no lower bound
    static const double NO_LOWER_BOUND = -2e19;

    /// @brief function type for functional f(x)
    typedef std::function<double(const Eigen::VectorXd&)> callback_f;

    /// Defines the optimization problems of the form
    ///  minₓ     f(x)       x ∈ Rⁿ
    ///           g(x) >= 0
    ///
    class OptimizationProblem {

    public:
        OptimizationProblem(const std::string& name);
        virtual ~OptimizationProblem();

        inline const std::string& name() const { return name_; }

        ////////////////////////////////////////////////////////////////////////
        // Objective function and its derivatives.

        /// @brief Evaulate the objective function.
        virtual double eval_f(const Eigen::VectorXd& x) = 0;

        /// @brief Evaulate the gradient of the objective function.
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) = 0;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        virtual Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x)
            = 0;

        ////////////////////////////////////////////////////////////////////////
        // Constraint function and its derivatives.

        /// @brief eval_g evaluates constraints at point x
        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) = 0;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        virtual Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) = 0;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        virtual std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x)
            = 0;

        ////////////////////////////////////////////////////////////////////////

        /// @brief Evaluate the objective and its derivatives.
        virtual void eval_f_and_fdiff(const Eigen::VectorXd& /*x*/,
            double& /*f_uk*/,
            Eigen::VectorXd& /*f_uk_jacobian*/)
        {
            throw NotImplementedError(
                "eval_f_and_fdiff level 1 not implemented");
        }

        /// @brief Evaluate the objective and its derivatives.
        virtual void eval_f_and_fdiff(const Eigen::VectorXd& /*x*/,
            double& /*f_uk*/,
            Eigen::VectorXd& /*f_uk_jacobian*/,
            Eigen::SparseMatrix<double>& /*f_uk_hessian*/)
        {
            throw NotImplementedError(
                "eval_f_and_fdiff level 2 not implemented");
        }

        virtual void eval_jac_g(const Eigen::VectorXd& /*x*/,
            Eigen::SparseMatrix<double>& /*jac_gx*/)
        {
            throw NotImplementedError("eval_jac_g sparse not implemented");
        }

        /// @brief eval_g evaluates constraints and jacobian at point x. Also
        /// returns list of active constraints (indices)
        virtual void eval_g(const Eigen::VectorXd& /*x*/,
            Eigen::VectorXd& /*gx*/,
            Eigen::SparseMatrix<double>& /*gx_jacobian*/,
            Eigen::VectorXi& /*gx_active*/)
        {
            throw NotImplementedError("eval_g level 3 not implemented");
        }

        /// @brief eval_g_and_gdiff evaluates constraints, jacobian and
        /// hessian at point x
        virtual void eval_g_and_gdiff(const Eigen::VectorXd& /*x*/,
            Eigen::VectorXd& /*gx*/,
            Eigen::MatrixXd& /*gx_jacobian*/,
            std::vector<Eigen::SparseMatrix<double>>& /*gx_hessian*/)
        {
            throw NotImplementedError("eval_g_and_gdiff not implemented");
        }

        ///////////////////////////////////////////////////////////////////////
        /// FINITE DIFFERENCES
        ///
        /// @brief Evaluate the gradient of the objective at point x using
        /// finitite differences.
        virtual Eigen::VectorXd eval_grad_f_approx(const Eigen::VectorXd& x);

        /// @brief Evaluate the hessian of the objective at point x using
        /// finitite differences.
        virtual Eigen::MatrixXd eval_hessian_f_approx(const Eigen::VectorXd& x);

        /// @brief Evaluate the jacobian of the constraint at point x using
        /// finitite differences.
        virtual Eigen::MatrixXd eval_jac_g_approx(const Eigen::VectorXd& x);

        /// @brief Evaluate the hessian of the constraint at point x using
        /// finitite differences.
        std::vector<Eigen::MatrixXd> eval_hessian_g_approx(
            const Eigen::VectorXd& /*x*/)
        {
            throw NotImplementedError("eval_hessian_g_approx not implemented");
        }

        virtual bool compare_grad_f_approx(
            const Eigen::VectorXd& x, const Eigen::VectorXd& grad);
        virtual bool compare_hessian_f_approx(const Eigen::VectorXd& x,
            const Eigen::SparseMatrix<double>& hessian);
        virtual bool compare_jac_g_approx(
            const Eigen::VectorXd& x, const Eigen::MatrixXd& jac);
        virtual bool compare_jac_g_approx(const Eigen::VectorXd& x,
            const Eigen::MatrixXd& jac,
            double& diff_norm);

        ///////////////////////////////////////////////////////////////////////
        /// FUNCTION POINTERS
        ///
        /// @brief creates function pointer that calls eval_f
        callback_f func_f();

        ///////////////////////////////////////////////////////////////////////
        /// BARRIER SPECIFIC
        ///
        virtual bool has_barrier_constraint() { return false; }
        virtual double get_barrier_epsilon()
        {
            throw NotImplementedError("get_barrier_epsilon not implemented");
        }
        virtual void set_barrier_epsilon(const double)
        {
            throw NotImplementedError("set_barrier_epsilon not implemented");
        }
        /// @brief Enable the line search mode.
        virtual void enable_line_search_mode(const Eigen::VectorXd& /*max_x*/)
        {
            throw NotImplementedError("set_barrier_epsilon not implemented");
        }
        /// @brief Disable the line search mode.
        virtual void disable_line_search_mode()
        {
            throw NotImplementedError("set_barrier_epsilon not implemented");
        }

        ////////////////////////////////////////////////////////////////////////
        /// @brief Call the intermediate_callback function.
        virtual bool eval_intermediate_callback(const Eigen::VectorXd& /*x*/)
        {
            return true;
        }

        /// indices if entry x_i is fixed
        virtual const Eigen::VectorXb& is_dof_fixed()
        {
            throw NotImplementedError("is_dof_fixed not implemented");
        }

        ////////////////////////////////////////////////////////////////////////
        // Fields
        int num_vars;        ///< @brief Number of variables
        int num_constraints; ///< @brief Number of constraints
        Eigen::VectorXd x0;  ///< @brief Initial value of x

        ////////////////////////////////////////////////////////////////////////
    private:
        std::string name_;
    };

} // namespace opt
} // namespace ccd
