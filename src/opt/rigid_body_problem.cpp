#include "rigid_body_problem.hpp"

#include <autodiff/finitediff.hpp>
#include <utils/flatten.hpp>
#include <utils/tensor.hpp>

namespace ccd {

namespace opt {
    RigidBodyProblem::RigidBodyProblem() {}
    RigidBodyProblem::~RigidBodyProblem() {}

    void RigidBodyProblem::initialize(
        ccd::physics::RigidBodySystem& rbs, CollisionConstraint& cstr)
    {
        model = std::make_shared<ccd ::physics::RigidBodySystem>(rbs);
        model->assemble();

        ParticlesDisplProblem::initialize(
            model->vertices, model->edges, model->displacements, cstr);
    }

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

    Eigen::VectorXd RigidBodyProblem::eval_grad_f_approx(
        const Eigen::VectorXd& x)
    {
        auto f = [&](const Eigen::VectorXd& v) -> double { return eval_f(v); };

        Eigen::VectorXd grad;
        ccd::finite_gradient(x, f, grad);
        return grad;
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

    Eigen::MatrixXd RigidBodyProblem::eval_hessian_f_approx(
        const Eigen::VectorXd& x)
    {

        auto f = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            return eval_grad_f(v);
        };

        Eigen::MatrixXd hess;
        ccd::finite_jacobian(x, f, hess);
        return hess;
    }
} // namespace opt
} // namespace ccd
