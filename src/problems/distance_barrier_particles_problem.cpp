#include "distance_barrier_particles_problem.hpp"

#include <autodiff/finitediff.hpp>
#include <utils/flatten.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

namespace opt {

    DistanceBarrierParticleProblem::DistanceBarrierParticleProblem(
        const std::string& name)
        : ParticlesProblem(name)
    {
    }

    void DistanceBarrierParticleProblem::settings(const nlohmann::json& params)
    {
        constraint_.settings(params["distance_barrier_constraint"]);
        opt_solver_.settings(params["barrier_solver"]);
        opt_solver_.set_problem(*this);
        ParticlesProblem::settings(params["particles_problem"]);
    }

    void DistanceBarrierParticleProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& xk,
        double& f_xk,
        Eigen::VectorXd& f_xk_grad,
        Eigen::SparseMatrix<double>& f_xk_hessian)
    {
        f_xk = eval_f(xk);
        f_xk_grad = eval_grad_f(xk);
        f_xk_hessian = eval_hessian_f(xk);
    }

    void DistanceBarrierParticleProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& sigma, double& f_uk, Eigen::VectorXd& f_uk_grad)
    {
        f_uk = eval_f(sigma);
        f_uk_grad = eval_grad_f(sigma);
    }

    void DistanceBarrierParticleProblem::update_constraints(
        const Eigen::MatrixXd& uk, const CstrSetFlag flag)
    {
        if (CstrSetFlag::UPDATE_CSTR_SET == flag) {
            constraint_.update_collision_set(uk);
        }
        constraint_.update_active_set(uk);
    }

    Eigen::VectorXd DistanceBarrierParticleProblem::eval_g(
        const Eigen::VectorXd& sigma)
    {
        return eval_g_set(sigma, CstrSetFlag::UPDATE_CSTR_SET);
    }

    Eigen::VectorXd DistanceBarrierParticleProblem::eval_g_set(
        const Eigen::VectorXd& xk, const CstrSetFlag flag)
    {
        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        update_constraints(uk, flag);

        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(uk, g_uk);
        return g_uk;
    }

    Eigen::MatrixXd DistanceBarrierParticleProblem::eval_jac_g(
        const Eigen::VectorXd& xk)
    {
        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        update_constraints(uk, CstrSetFlag::UPDATE_CSTR_SET);

        Eigen::MatrixXd g_uk_jacobian;
        constraint_.compute_constraints_jacobian(uk, g_uk_jacobian);

#ifdef WITH_DERIVATIVE_CHECK
        compare_jac_g_approx(xk, g_uk_jacobian);
#endif
        return g_uk_jacobian;
    }

    std::vector<Eigen::SparseMatrix<double>>
    DistanceBarrierParticleProblem::eval_hessian_g(const Eigen::VectorXd& xk)
    {
        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        update_constraints(uk, CstrSetFlag::UPDATE_CSTR_SET);

        std::vector<Eigen::SparseMatrix<double>> g_uk_hessian;
        constraint_.compute_constraints_hessian(uk, g_uk_hessian);
        return g_uk_hessian;
    }

    void DistanceBarrierParticleProblem::eval_g_and_gdiff(
        const Eigen::VectorXd& xk,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian)
    {
        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        update_constraints(uk, CstrSetFlag::UPDATE_CSTR_SET);

        constraint_.compute_constraints(uk, g_uk);
        constraint_.compute_constraints_jacobian(uk, g_uk_jacobian);
        constraint_.compute_constraints_hessian(uk, g_uk_hessian);

#ifdef WITH_DERIVATIVE_CHECK
        compare_jac_g_approx(xk, g_uk_jacobian);
#endif
    }

    void DistanceBarrierParticleProblem::compare_jac_g_approx(
        const Eigen::MatrixXd& xk, const Eigen::MatrixXd& g_uk_jacobian)
    {

        Eigen::MatrixXd jac_g_approx = eval_jac_g_approx(*this, xk);
        if (!compare_jacobian(g_uk_jacobian, jac_g_approx)) {
            spdlog::error(
                "barrier_particles status=fail "
                "message='constraint jacobian finite-differences failed'");
        }
    }

} // namespace opt
} // namespace ccd
