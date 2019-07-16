#include "ad_hoc_problem.hpp"

#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace opt {

    AdHocProblem::AdHocProblem()
        : OptimizationProblem("AdHoc")
    {
        this->f = [](const Eigen::VectorXd&) -> double {
            throw NotImplementedError("Objective function not implemented!");
        };
        this->grad_f = [](const Eigen::VectorXd&) -> Eigen::VectorXd {
            throw NotImplementedError(
                "Gradient of the objective function not implemented!");
        };
        this->hessian_f = [](const Eigen::VectorXd&) -> Eigen::MatrixXd {
            throw NotImplementedError(
                "Hessian of the objective function not implemented!");
        };
        this->g = [](const Eigen::VectorXd&) -> Eigen::VectorXd {
            throw NotImplementedError("Constraint function not implemented!");
        };
        this->jac_g = [](const Eigen::VectorXd&) -> Eigen::MatrixXd {
            throw NotImplementedError(
                "Jacobian of the constraint function not implemented!");
        };
        this->hessian_g = [](const Eigen::VectorXd&)
            -> std::vector<Eigen::SparseMatrix<double>> {
            throw NotImplementedError("Second derivative of the constraint "
                                      "function not implemented!");
        };
    }

    AdHocProblem::AdHocProblem(int num_vars, int num_constraints)
        : AdHocProblem()
    {
        this->num_vars = num_vars;
        this->num_constraints = num_constraints;

        this->x0.resize(num_vars);
        this->x0.setConstant(0.0);
    }

    double AdHocProblem::eval_f(const Eigen::VectorXd& x) { return f(x); }

    Eigen::VectorXd AdHocProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return grad_f(x);
    }
    Eigen::SparseMatrix<double> AdHocProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        return hessian_f(x).sparseView();
    }


    void AdHocProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& x, double& value, Eigen::VectorXd& grad)
    {
        value = eval_f(x);
        grad = eval_grad_f(x);
    }

    void AdHocProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& value,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hessian)
    {
        value = eval_f(x);
        grad = eval_grad_f(x);
        hessian = eval_hessian_f(x);
    }


    Eigen::VectorXd AdHocProblem::eval_g(const Eigen::VectorXd& x)
    {
        return g(x);
    };

    Eigen::MatrixXd AdHocProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        return jac_g(x);
    };

    std::vector<Eigen::SparseMatrix<double>> AdHocProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        return hessian_g(x);
    };

    void AdHocProblem::eval_jac_g(
        const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& jac_gx)
    {
        jac_gx = jac_g(x).sparseView();
    };

    void AdHocProblem::eval_g(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::SparseMatrix<double>& gx_jacobian,
        Eigen::VectorXi& gx_active)
    {
        gx = eval_g(x);
        eval_jac_g(x, gx_jacobian);
        gx_active = Eigen::VectorXi::LinSpaced(gx.rows(), 0, int(gx.rows()));
    }

    void AdHocProblem::eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& gx_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        gx = eval_g(x);
        gx_jacobian = eval_jac_g(x);
        gx_hessian = eval_hessian_g(x);
    }

} // namespace opt
} // namespace ccd
