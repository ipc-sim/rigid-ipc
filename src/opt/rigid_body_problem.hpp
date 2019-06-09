#pragma once

#include <opt/particles_problem.hpp>
#include <physics/rigid_body_system.hpp>

namespace ccd {

namespace opt {

    class RigidBodyProblem : public ParticlesDisplProblem {
    public:
        RigidBodyProblem();
        ~RigidBodyProblem() override;

        void initialize(ccd::physics::RigidBodySystem& rbs,
            CollisionConstraint& cstr);

        double eval_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override;
        Eigen::VectorXd eval_grad_f_approx(const Eigen::VectorXd& x);

        Eigen::SparseMatrix<double> eval_hessian_f_sparse(
            const Eigen::VectorXd& x) override;
        Eigen::MatrixXd eval_hessian_f_approx(
            const Eigen::VectorXd& x);
        std::shared_ptr<ccd::physics::RigidBodySystem> model;

    };
} // namespace opt
} // namespace ccd
