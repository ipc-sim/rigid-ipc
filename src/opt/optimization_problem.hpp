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

    class IUnconstraintedProblem {
    public:
        virtual ~IUnconstraintedProblem() = default;

        /// @brief Evaulate the objective function.
        virtual double eval_f(const Eigen::VectorXd& x) = 0;

        /// @brief Evaulate the gradient of the objective function.
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) = 0;

        /// @brief Evaluate the hessian of the objective as a sparse matrix.
        virtual Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x)
            = 0;

        virtual const int& num_vars() = 0;
        virtual const Eigen::VectorXd& starting_point() = 0;
    };

    /// Defines  optimization problems of the form
    ///  minₓ     f(x)       x ∈ Rⁿ
    ///           g(x) >= 0
    ///

    enum class CstrSetFlag { UPDATE_CSTR_SET = 0, KEEP_CSTR_SET };

    class IConstraintedProblem : public virtual IUnconstraintedProblem {
    public:
        virtual ~IConstraintedProblem() = default;

        /// @brief eval_g evaluates constraints at point x
        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) = 0;

        /// \brief eval_g_set: evals the constraints with option of updating or
        /// freezing constraint set
        ///
        virtual Eigen::VectorXd eval_g_set(
            const Eigen::VectorXd& x, const CstrSetFlag flag)
            = 0;

        /// @brief eval_jac_g evaluates constraints jacobian at point x
        virtual Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) = 0;

        // @brief eval_hessian_g evaluates constraints hessian at point x
        virtual std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x)
            = 0;

        virtual const int& num_constraints() = 0;
    };

    //  Used for the Volume Constraint Solvers (i.e just NCP)
    class INCPProblem : public virtual IUnconstraintedProblem {
    public:
        virtual ~INCPProblem() = default;

        /// @brief eval_g evaluates constraints at point x
        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) = 0;
        virtual void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active)
            = 0;

        virtual const Eigen::VectorXb& is_dof_fixed() = 0;

    };

    class IBarrierGeneralProblem : public virtual IConstraintedProblem {
    public:
        virtual ~IBarrierGeneralProblem() = default;

        virtual void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian)
            = 0;
        virtual void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian)
            = 0;

        virtual void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
            = 0;

        virtual double get_barrier_epsilon() = 0;
        virtual void set_barrier_epsilon(const double eps) = 0;
        virtual const Eigen::VectorXb& is_dof_fixed() = 0;
    };

    class IBarrierProblem : public virtual IUnconstraintedProblem {
    public:
        virtual ~IBarrierProblem() = default;

        virtual double eval_f_set(
            const Eigen::VectorXd& x, const CstrSetFlag flag)
            = 0;

        virtual void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian)
            = 0;

        virtual void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian)
            = 0;

        virtual double get_barrier_epsilon() = 0;
        virtual const Eigen::VectorXb& is_dof_fixed() = 0;
    };

    /// Helper Functions for checking finite differences
    /// ------------------------------------------------
    Eigen::VectorXd eval_grad_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x);
    Eigen::MatrixXd eval_hess_f_approx(
        IUnconstraintedProblem& problem, const Eigen::VectorXd& x);

    Eigen::MatrixXd eval_jac_g_approx(
        IConstraintedProblem& problem, const Eigen::VectorXd& x);

} // namespace opt
} // namespace ccd
