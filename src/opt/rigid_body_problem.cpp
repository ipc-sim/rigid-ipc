#include "rigid_body_problem.hpp"

#include <utils/flatten.hpp>
#include <utils/tensor.hpp>

namespace ccd {

namespace opt {
    RigidBodyProblem::RigidBodyProblem()
        : model(nullptr)
    {
    }

    RigidBodyProblem::~RigidBodyProblem() {}

    void RigidBodyProblem::initialize(
        ccd::physics::RigidBodySystem& rbs, CollisionConstraint& cstr)
    {
        model = std::make_shared<ccd ::physics::RigidBodySystem>(rbs);
        model->assemble();

        init_num_vars();
        fixed_dof = Eigen::MatrixXb::Zero(num_vars, 1);

        ParticlesDisplProblem::initialize(
            model->vertices, model->edges, model->displacements, cstr);
    }

    void RigidBodyProblem::init_num_vars()
    {
        num_vars = 3 * this->model->get_rigid_body_count();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Objective function and its derivatives.

    double RigidBodyProblem::eval_f(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);
        return ParticlesDisplProblem::eval_f(Uk);
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        // TODO: compute displacement, gradient and hessian together.
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);

        Eigen::SparseMatrix<double> jac_uk_x;
        model->compute_displacements_gradient(x, jac_uk_x);

        Eigen::VectorXd grad_uk = ParticlesDisplProblem::eval_grad_f(Uk);
        return grad_uk.transpose() * jac_uk_x;
    }

    Eigen::MatrixXd RigidBodyProblem::eval_hessian_f(const Eigen::VectorXd& x)
    {
        return Eigen::MatrixXd(eval_hessian_f_sparse(x));
    }

    Eigen::SparseMatrix<double> RigidBodyProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        // TODO: compute displacement, gradient and hessian together.
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);

        std::vector<Eigen::SparseMatrix<double>> hess_uk_x;
        model->compute_displacements_hessian(x, hess_uk_x);

        Eigen::SparseMatrix<double> jac_uk_x;
        model->compute_displacements_gradient(x, jac_uk_x);

        Eigen::SparseMatrix<double> hess_uk
            = ParticlesDisplProblem::eval_hessian_f_sparse(Uk);
        Eigen::VectorXd grad_uk = ParticlesDisplProblem::eval_grad_f(Uk);

        // chain rule
        return jac_uk_x.transpose() * hess_uk * jac_uk_x
            + tensor::multiply(grad_uk, hess_uk_x);
    }
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Constraint function and its derivatives.

    Eigen::VectorXd RigidBodyProblem::eval_g(const Eigen::VectorXd& x)
    {
        // Compute the displacements.
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);
        // Compute the gradient of the displacements g(U(x))
        return ParticlesDisplProblem::eval_g(displacements);
    }

    // Evaluate the constraint's jacobian at point x.
    Eigen::MatrixXd RigidBodyProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        // Compute the displacements, U.
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        // Compute ∇U(x)
        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        // Compute ∇g(U)
        Eigen::MatrixXd g_gradient
            = ParticlesDisplProblem::eval_jac_g(displacements);

        // Apply the chain rule to compute ∇g(U(x)) = ∇g(U) * ∇U(x)
        return g_gradient * displacements_gradient;
    }

    // Evaluate the constraint's hessian at point x.
    std::vector<Eigen::SparseMatrix<double>> RigidBodyProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        // Compute the displacements U
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        // Compute ∇U(x)
        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        // Compute ∇²U(x)
        std::vector<Eigen::SparseMatrix<double>> displacements_hessian;
        model->compute_displacements_hessian(x, displacements_hessian);

        Eigen::VectorXd _;                                  // g(U(x))
        Eigen::MatrixXd g_gradient;                         // ∇g(U)
        std::vector<Eigen::SparseMatrix<double>> g_hessian; // ∇²g(U)
        ParticlesDisplProblem::eval_g_and_gdiff(
            displacements, _, g_gradient, g_hessian);

        // Compute ∇²g(U(x))
        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        gx_hessian.reserve(g_hessian.size());
        for (unsigned long i = 0; i < g_hessian.size(); i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(displacements_gradient.transpose()
                    * g_hessian[i] * displacements_gradient
                + tensor::multiply(
                    g_gradient.row(long(i)), displacements_hessian));
        }

        return gx_hessian;
    }

    // Evaluates the constraints, jacobian, and hessian at point x.
    void RigidBodyProblem::eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx, Eigen::MatrixXd& gx_gradient,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        std::vector<Eigen::SparseMatrix<double>> displacements_hessian;
        model->compute_displacements_hessian(x, displacements_hessian);

        Eigen::MatrixXd g_gradient;
        std::vector<Eigen::SparseMatrix<double>> g_hessian;
        ParticlesDisplProblem::eval_g_and_gdiff(
            displacements, gx, g_gradient, g_hessian);

        gx_gradient = g_gradient * displacements_gradient;

        gx_hessian.clear();
        gx_hessian.reserve(g_hessian.size());
        for (unsigned long i = 0; i < g_hessian.size(); i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(displacements_gradient.transpose()
                    * g_hessian[i] * displacements_gradient
                + tensor::multiply(
                    g_gradient.row(long(i)), displacements_hessian));
        }
    }
    ////////////////////////////////////////////////////////////////////////////

} // namespace opt
} // namespace ccd
