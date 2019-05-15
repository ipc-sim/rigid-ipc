#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace opt {

    ///@brief default value for no upper bound
    static const double NO_UPPER_BOUND = 2e19;
    ///@brief default value for no lower bound
    static const double NO_LOWER_BOUND = -2e19;

    ///@brief function type for functional f(x)
    typedef std::function<double(const Eigen::VectorXd&)> callback_f;
    ///@brief function type for gradient of functional \nabla f(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd&)>
        callback_grad_f;
    ///@brief function type for hessian of functional \f$\nabla^2 f(x)\f$
    typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>
        callback_hessian_f;
    ///@brief function type for constraints g(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd&)> callback_g;
    ///@brief function type for jacobian of constraints \nabla g(x)
    typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>
        callback_jac_g;
    ///@brief function type for direvative of the jacobian of constraints
    /// \f$\nabla^2 g(x)\f$
    typedef std::function<std::vector<Eigen::SparseMatrix<double>>(
        const Eigen::VectorXd&)>
        callback_hessian_g;

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
    /**
     *  Defines the optimization problems of the form
     *      MIN     f(x)      x ∈ Rⁿ
     *
     *   s.t.       g_L ≤ g(x) ≤ g_U
     *              x_L ≤  x   ≤ x_U
     */

    class OptimizationProblem {
    public:
        virtual ~OptimizationProblem();

        virtual double eval_f(const Eigen::VectorXd& x) = 0;
        virtual Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) = 0;
        virtual Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) = 0;
        virtual Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& x);

        virtual Eigen::VectorXd eval_g(const Eigen::VectorXd& x) = 0;
        virtual Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) = 0;
        virtual std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x)
            = 0;

        callback_f func_f();
        callback_grad_f func_grad_f();
        callback_g func_g();

        /// @brief Check that the problem is valid and initalized
        bool validate_problem();

        /// @brief Check if all constraints are satisfied at a location.
        bool are_constraints_satisfied(
            const Eigen::VectorXd& x, const double tol);

        /// @brief fixed_dof[i] == true indicates x[i] is not a variable
        ArrayXb fixed_dof;

        int num_vars;            ///< @brief Number of variables
        int num_constraints;     ///< @brief Number of constraints
        Eigen::VectorXd x0;      ///< @brief Initial value of x
        Eigen::VectorXd x_lower; ///< @brief Lower bound of x
        Eigen::VectorXd x_upper; ///< @brief Upper bound of x
        Eigen::VectorXd g_lower; ///< @brief Lower bound of the constraint
        Eigen::VectorXd g_upper; ///< @brief Upper bound of the constraint
    };

    class AdHocProblem : public OptimizationProblem {
    public:
        callback_f f;           ///< @brief Objective function
        callback_grad_f grad_f; ///< @brief Gradient of the objective function
        callback_hessian_f
            hessian_f;        ///< @brief Hessian of the objective function
        callback_g g;         ///< @brief Constraint function
        callback_jac_g jac_g; ///< @brief Jacobian of the constraint function
        callback_hessian_g
            hessian_g; ///< @brief Hessian of the constraint function

        /// @brief Default constructor
        AdHocProblem();

        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        /// @brief Resize fields accordingly
        AdHocProblem(int num_vars, int num_constraints);
    };

} // namespace opt
} // namespace ccd