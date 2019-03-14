#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ccd {
namespace opt {

    /**
     *  We handle optimization problems of the form
     *      MIN     f(x)      x ∈ Rⁿ
     *
     *   s.t.       g_L ≤ g(x) ≤ g_U
     *              x_L ≤  x   ≤ x_U
     *
     */

    ///@brief default value for no upper bound
    static const double NO_UPPER_BOUND = 2e19;
    ///@brief default value for no lower bound
    static const double NO_LOWER_BOUND = -2e19;

    ///@brief function type for functional f(x)
    typedef std::function<double(const Eigen::VectorXd& x)> callback_f;
    ///@brief function type for gradient of functional \nabla f(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>
        callback_grad_f;
    ///@brief function type for constraints g(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd& x)> callback_g;
    ///@brief function type for jacobian of constraints \nabla g(x)
    typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd& x)>
        callback_jac_g;
    ///@brief function type for callback, called once on each iteration
    typedef std::function<void(const Eigen::VectorXd& x, const double obj_value,
        const Eigen::VectorXd& dual, const int iteration)>
        callback_intermediate;

    /// @brief Methods available for
    enum OptimizationMethod {
        MMA,   ///<@brief Method of Moving Asymptotes (NLopt)
        SLSQP, ///<@brief Sequential Least-Squares Quadratic Programming (NLopt)
        IP     ///<@brief Interior-Point Method (Ipopt)
    };

    struct OptimizationResult {
        Eigen::MatrixXd x; ///<@brief the solution of the optimization.
        double fun;        ///<@brief value of the objective function
        bool success;      ///<@brief whether or not the optimizer exited
                           ///< successfully.
    };

    struct OptimizationProblem {
        int num_vars;
        int num_constraints;
        Eigen::VectorXd x0;
        Eigen::VectorXd x_lower;
        Eigen::VectorXd x_upper;
        Eigen::VectorXd g_lower;
        Eigen::VectorXd g_upper;
        callback_f f;
        callback_grad_f grad_f;
        callback_g g;
        callback_jac_g jac_g;
        int verbosity = 0;
        int max_iter = 3000;
        callback_intermediate intermediate_cb;

        OptimizationProblem();
        OptimizationProblem(const Eigen::VectorXd& x0, const callback_f f,
            const callback_grad_f grad_f,
            const Eigen::VectorXd& x_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& x_upper = Eigen::VectorXd(),
            const int num_constraints = 0, const callback_g& g = nullptr,
            const callback_jac_g& jac_g = nullptr,
            const Eigen::VectorXd& g_lower = Eigen::VectorXd(),
            const Eigen::VectorXd& g_upper = Eigen::VectorXd(),
            const callback_intermediate callback = nullptr);
        bool validate_problem();
    };

    OptimizationResult minimize(
        const OptimizationMethod& method, OptimizationProblem& problem);

} // namespace opt

} // namespace ccd
