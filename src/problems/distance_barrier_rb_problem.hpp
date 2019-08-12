#pragma once

#include <autodiff/autodiff_types.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <physics/rigid_body_problem.hpp>

namespace ccd {

class DistanceBarrierRBProblem : public physics::RigidBodyProblem,
                                 public virtual opt::IBarrierGeneralProblem {
public:
    ////////////////////////////////////////////////////////////////
    /// IConstrainedProblem

    /// @brief eval_g evaluates constraints at point x
    Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override;

    /// @brief eval_jac_g evaluates constraints jacobian at point x
    Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override;

    // @brief eval_hessian_g evaluates constraints hessian at point x
    std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
        const Eigen::VectorXd& x) override;

    ////////////////////////////////////////////////////////////
    /// IBarrierProblem

    Eigen::VectorXd eval_g_(
        const Eigen::VectorXd& x, const bool update_constraint_set) override;

    void eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk,
        Eigen::VectorXd& f_uk_jacobian,
        Eigen::SparseMatrix<double>& f_uk_hessian) override;

    void eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk,
        Eigen::VectorXd& f_uk_jacobian) override;

    void eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& gx_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian) override;

    double get_barrier_epsilon() override;
    void set_barrier_epsilon(const double eps) override;
    const Eigen::VectorXb& is_dof_fixed() override;

protected:
    void update_constraints(
        const Eigen::MatrixXd& xk, const bool update_cstr_set);

    opt::DistanceBarrierConstraint constraint_;
};

} // namespace ccd
