#include <opt/optimization_problem.hpp>

namespace ccd {
namespace opt {

    /// @brief function type for gradient of functional \nabla f(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd&)>
        callback_grad_f;
    /// @brief function type for hessian of functional \f$\nabla^2 f(x)\f$
    typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>
        callback_hessian_f;
    /// @brief function type for constraints g(x)
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd&)> callback_g;
    /// @brief function type for jacobian of constraints \nabla g(x)
    typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>
        callback_jac_g;
    /// @brief function type for direvative of the jacobian of constraints
    /// \f$\nabla^2 g(x)\f$
    typedef std::function<std::vector<Eigen::SparseMatrix<double>>(
        const Eigen::VectorXd&)>
        callback_hessian_g;

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
        std::function<Eigen::VectorXb&()> fis_dof_fixed;

        /// @brief Default constructor
        AdHocProblem();
        /// @brief Resize fields accordingly
        AdHocProblem(int num_vars, int num_constraints);

        /// FUNCTIONAL
        /// /////////////////////////
        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override;

        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& value,
            Eigen::VectorXd& grad) override;
        void eval_f_and_fdiff(const Eigen::VectorXd& /*x*/,
            double& /*f_uk*/,
            Eigen::VectorXd& /*f_uk_jacobian*/,
            Eigen::SparseMatrix<double>& /*f_uk_hessian*/) override;

        /// /////////////////////////
        /// CONSTRAINTS
        /// /////////////////////////
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override;

        virtual void eval_jac_g(const Eigen::VectorXd& x,
            Eigen::SparseMatrix<double>& jac_gx) override;
        virtual void eval_g(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::SparseMatrix<double>& gx_jacobian,
            Eigen::VectorXi& gx_active) override;

        virtual void eval_g_and_gdiff(const Eigen::VectorXd& /*x*/,
            Eigen::VectorXd& /*gx*/,
            Eigen::MatrixXd& /*gx_jacobian*/,
            std::vector<Eigen::SparseMatrix<double>>& /*gx_hessian*/) override;

        const Eigen::VectorXb& is_dof_fixed() override
        {
            return fis_dof_fixed();
        }

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
