#include <opt/optimization_problem.hpp>

namespace ccd {
namespace opt {

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

        std::function<bool()> fhas_barrier_constraint;
        std::function<double()> fget_barrier_epsilon;
        std::function<void(const double)> fset_barrier_epsilon;

        /// @brief Default constructor
        AdHocProblem();
        /// @brief Resize fields accordingly
        AdHocProblem(int num_vars, int num_constraints);

        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_hessian_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        bool has_barrier_constraint() override
        {
            if (fhas_barrier_constraint != nullptr) {
                return fhas_barrier_constraint();
            }
            return OptimizationProblem::has_barrier_constraint();
        }

        double get_barrier_epsilon() override
        {
            if (fget_barrier_epsilon != nullptr) {
                return fget_barrier_epsilon();
            }
            return OptimizationProblem::get_barrier_epsilon();
        }
        void set_barrier_epsilon(const double eps) override
        {
            if (fset_barrier_epsilon != nullptr) {
                return fset_barrier_epsilon(eps);
            }
            return OptimizationProblem::set_barrier_epsilon(eps);
        }
    };

} // namespace opt
} // namespace ccd
