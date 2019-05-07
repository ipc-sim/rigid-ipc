#pragma once

#include <Eigen/Core>
#include <vector>

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
    typedef std::function<std::vector<Eigen::MatrixXd>(const Eigen::VectorXd&)>
        callback_hessian_g;

    /**
     *  Defines the optimization problems of the form
     *      MIN     f(x)      x ∈ Rⁿ
     *
     *   s.t.       g_L ≤ g(x) ≤ g_U
     *              x_L ≤  x   ≤ x_U
     *
     */
    struct OptimizationProblem {
        int num_vars;            ///< @brief Number of variables
        int num_constraints;     ///< @brief Number of constraints
        Eigen::VectorXd x0;      ///< @brief Initial value of x
        Eigen::VectorXd x_lower; ///< @brief Lower bound of x
        Eigen::VectorXd x_upper; ///< @brief Upper bound of x
        Eigen::VectorXd g_lower; ///< @brief Lower bound of the constraint
        Eigen::VectorXd g_upper; ///< @brief Upper bound of the constraint
        callback_f f;            ///< @brief Objective function
        callback_grad_f grad_f;  ///< @brief Gradient of the objective function
        callback_hessian_f
            hessian_f;        ///< @brief Hessian of the objective function
        callback_g g;         ///< @brief Constraint function
        callback_jac_g jac_g; ///< @brief Jacobian of the constraint function
        callback_hessian_g
            hessian_g; ///< @brief Hessian of the constraint function
        Eigen::Array<bool, Eigen::Dynamic, 1> fixed_dof;

        /// @brief Default constructor
        OptimizationProblem();
        /// @brief Resize fields accordingly
        OptimizationProblem(int num_vars, int num_constraints);
        /// @brief Construct an optimization problem
        OptimizationProblem(const Eigen::VectorXd& x0, const callback_f f,
            const callback_grad_f grad_f,
            const Eigen::VectorXd& x_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& x_upper = Eigen::VectorXd(),
            const int num_constraints = 0, const callback_g& g = nullptr,
            const callback_jac_g& jac_g = nullptr,
            const Eigen::VectorXd& g_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& g_upper = Eigen::VectorXd());
        /// @brief Construct an optimization problem with hessians
        OptimizationProblem(const Eigen::VectorXd& x0, const callback_f f,
            const callback_grad_f grad_f,
            const callback_hessian_f& hessian_f = nullptr,
            const Eigen::VectorXd& x_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& x_upper = Eigen::VectorXd(),
            const int num_constraints = 0, const callback_g& g = nullptr,
            const callback_jac_g& jac_g = nullptr,
            const callback_hessian_g& hessian_g = nullptr,
            const Eigen::VectorXd& g_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& g_upper = Eigen::VectorXd());
        /// @brief Check that the problem is valid and initalized
        bool validate_problem();
        /// @brief Check if all constraints are satisfied at a location.
        bool are_constraints_satisfied(
            const Eigen::VectorXd& x, const double tol) const;
    };

} // namespace opt
} // namespace ccd
