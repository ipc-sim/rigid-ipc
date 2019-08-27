#include "distance_barrier_rb_problem.hpp"

#include <iostream>
#include <utils/tensor.hpp>

#include <autodiff/finitediff.hpp>
#include <logger.hpp>
#include <profiler.hpp>

#include <complex.h>
#include <tgmath.h>

namespace ccd {

namespace opt {

    DistanceBarrierRBProblem::DistanceBarrierRBProblem(const std::string& name)
        : RigidBodyProblem(name)
    {
    }

    void DistanceBarrierRBProblem::settings(const nlohmann::json& params)
    {
        constraint_.settings(params["distance_barrier_constraint"]);
        opt_solver_.settings(params["barrier_solver"]);
        opt_solver_.set_problem(*this);
        RigidBodyProblem::settings(params["rigid_body_problem"]);
    }

    void DistanceBarrierRBProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& sigma,
        double& f_uk,
        Eigen::VectorXd& f_uk_grad,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        f_uk = eval_f(sigma);
        f_uk_grad = eval_grad_f(sigma);
        f_uk_hessian = eval_hessian_f(sigma);
    }

    void DistanceBarrierRBProblem::eval_f_and_fdiff(
        const Eigen::VectorXd& sigma, double& f_uk, Eigen::VectorXd& f_uk_grad)
    {
        f_uk = eval_f(sigma);
        f_uk_grad = eval_grad_f(sigma);
    }

    Eigen::VectorXd DistanceBarrierRBProblem::eval_g(
        const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(uk, g_uk);
        return g_uk;
    }

    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g(
        const Eigen::VectorXd& sigma)
    {
        NAMED_PROFILE_POINT("eval_jac_g__update_constraints", UPDATE)
        NAMED_PROFILE_POINT("eval_jac_g__eval_jac", EVAL)

        PROFILE_START(UPDATE)

        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        EdgeVertexCandidates ev_candidates;
        auto check = constraint_.get_active_barrier_set(uk, ev_candidates);
        assert(check == DistanceBarrierConstraint::NO_COLLISIONS);
        PROFILE_END(UPDATE)

        PROFILE_START(EVAL)
        Eigen::MatrixXd gx_jacobian = eval_jac_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL)

        bool derivative_check = true;
        if (derivative_check) {
            compare_jac_g(sigma, ev_candidates, gx_jacobian);
        }
        return gx_jacobian;
    }
    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g_full(
        const Eigen::VectorXd& sigma, const EdgeVertexCandidates& ev_candidates)
    {

        DiffScalarBase::setVariableCount(size_t(num_vars_));

        // create the variables
        typedef DScalar1<double, Eigen::VectorXd> DDouble1;
        typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, 1> D1VectorXd;
        typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, Eigen::Dynamic>
            D1MatrixXd;
        D1VectorXd d_sigma;
        d_sigma.resize(sigma.rows());
        for (int r = 0; r < sigma.rows(); r++) {
            d_sigma[r] = DDouble1(size_t(r), sigma[r]);
        }

        D1VectorXd d_qk
            = m_assembler.m_dof_to_position.cast<DDouble1>() * d_sigma;
        D1MatrixXd d_uk = m_assembler.world_vertices<DDouble1>(d_qk)
            - vertices_t0.cast<DDouble1>();
        D1VectorXd d_g_uk;
        constraint_.compute_candidates_constraints<DDouble1>(
            d_uk, ev_candidates, d_g_uk);

        Eigen::MatrixXd jac_g;
        jac_g.resize(ev_candidates.size(), num_vars_);
        for (int i = 0; i < d_g_uk.rows(); i++) {
            jac_g.row(i) = d_g_uk(i).getGradient();
        }

        return jac_g;
    }

    std::vector<Eigen::SparseMatrix<double>>
    DistanceBarrierRBProblem::eval_hessian_g(const Eigen::VectorXd& sigma)
    {
        NAMED_PROFILE_POINT("eval_hess_g__update_constraints", UPDATE)
        NAMED_PROFILE_POINT("eval_hess_g__eval", EVAL)

        PROFILE_START(UPDATE)
        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        EdgeVertexCandidates ev_candidates;
        auto check = constraint_.get_active_barrier_set(uk, ev_candidates);
        assert(check == DistanceBarrierConstraint::NO_COLLISIONS);

        PROFILE_END(UPDATE)

        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        PROFILE_START(EVAL)
        gx_hessian = eval_hessian_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL)

        return gx_hessian;
    }

    void DistanceBarrierRBProblem::eval_g_and_gdiff(
        const Eigen::VectorXd& sigma,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& gx_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        NAMED_PROFILE_POINT("eval_g_and_gdiff__update_constraints", UPDATE)
        NAMED_PROFILE_POINT("eval_hess_g__eval_grad", EVAL_GRAD)
        NAMED_PROFILE_POINT("eval_hess_g__eval_hess", EVAL_HESS)

        PROFILE_START(UPDATE)
        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        EdgeVertexCandidates ev_candidates;
        auto check = constraint_.get_active_barrier_set(uk, ev_candidates);
        assert(check == DistanceBarrierConstraint::NO_COLLISIONS);
        PROFILE_END(UPDATE)

        constraint_.compute_candidates_constraints(uk, ev_candidates, gx);

        PROFILE_START(EVAL_GRAD)
        gx_jacobian = eval_jac_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL_GRAD)

        PROFILE_START(EVAL_HESS)
        gx_hessian = eval_hessian_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL_HESS)

        bool derivative_check = true;
        if (derivative_check) {

            assert(compare_jac_g(sigma, ev_candidates, gx_jacobian));
        }
    }

    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g_core(
        const Eigen::VectorXd& sigma,
        const EdgeVertexCandidates& distance_candidates)
    {
        Eigen::MatrixXd jac_g;
        jac_g.resize(distance_candidates.size(), num_vars_);
        jac_g.setZero();

        typedef AutodiffType<6> Diff;

        for (size_t i = 0; i < distance_candidates.size(); ++i) {
            const auto& ev_candidate = distance_candidates[i];

            RB2Candidate rbc;
            extract_local_system(ev_candidate, rbc);

            Eigen::VectorXd gradient
                = distance_barrier<Diff::DDouble1>(sigma, rbc).getGradient();

            size_t cstr_id = i;
            jac_g.block(int(cstr_id), 3 * rbc.vertex_body_id, 1, 3)
                = gradient.segment(0, 3).transpose();
            jac_g.block(int(cstr_id), 3 * rbc.edge_body_id, 1, 3)
                = gradient.segment(3, 3).transpose();
        }

        return jac_g;
    }

    std::vector<Eigen::SparseMatrix<double>>
    DistanceBarrierRBProblem::eval_hessian_g_core(const Eigen::VectorXd& sigma,
        const EdgeVertexCandidates& distance_candidates)
    {
        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        gx_hessian.resize(distance_candidates.size());

        typedef AutodiffType<6> Diff;
        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        for (size_t i = 0; i < distance_candidates.size(); ++i) {
            const auto& ev_candidate = distance_candidates[i];

            RB2Candidate rbc;
            extract_local_system(ev_candidate, rbc);

            Eigen::MatrixXd hessian
                = distance_barrier<Diff::DDouble2>(sigma, rbc).getHessian();

            triplets.clear();
            triplets.reserve(6 * 6);

            int bodies[2] = { rbc.vertex_body_id, rbc.edge_body_id };

            for (int b_i = 0; b_i < 2; b_i++) {
                for (int b_j = 0; b_j < 2; b_j++) {
                    for (int dim_i = 0; dim_i < 3; dim_i++) {
                        for (int dim_j = 0; dim_j < 3; dim_j++) {
                            double v
                                = hessian(3 * b_i + dim_i, 3 * b_j + dim_j);
                            int r = 3 * bodies[b_i] + dim_i;
                            int c = 3 * bodies[b_j] + dim_j;
                            triplets.push_back(M(r, c, v));
                        }
                    }
                }
            }
            Eigen::SparseMatrix<double> global_el_hessian(num_vars_, num_vars_);
            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());

            size_t cstr_id = i;
            gx_hessian[cstr_id] = global_el_hessian;
        }

        return gx_hessian;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance_barrier(
        const Eigen::VectorXd& sigma, const RB2Candidate& rbc)
    {
        T d = distance<T>(sigma, rbc);
        T barrier = constraint_.distance_barrier<T>(d);
        return barrier;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RB2Candidate& rbc)
    {
        typedef AutodiffType<6> Diff;

        Diff::activate();
        typedef Eigen::Matrix<T, 3, 1> DTVector3d;
        typedef Eigen::Matrix<T, 2, 1> DTVector2d;

        DTVector3d sigma_E, sigma_V, position_E, position_V;

        sigma_V = Diff::dTvars<T>(0, sigma.segment(3 * rbc.vertex_body_id, 3));
        sigma_E = Diff::dTvars<T>(3, sigma.segment(3 * rbc.edge_body_id, 3));

        position_V = sigma_V.array()
            * m_assembler.m_dof_to_position.diagonal()
                  .segment(3 * rbc.vertex_body_id, 3)
                  .cast<T>()
                  .array();

        position_E = sigma_E.array()
            * m_assembler.m_dof_to_position.diagonal()
                  .segment(3 * rbc.edge_body_id, 3)
                  .cast<T>()
                  .array();

        const auto& rbs = m_assembler.m_rbs;
        DTVector2d da = rbs[size_t(rbc.edge_body_id)].world_vertex<T>(
            position_E, rbc.edge0_local_id);
        DTVector2d db = rbs[size_t(rbc.edge_body_id)].world_vertex<T>(
            position_E, rbc.edge1_local_id);
        DTVector2d dc = rbs[size_t(rbc.vertex_body_id)].world_vertex<T>(
            position_V, rbc.vertex_local_id);

        T distance = sqrt(point_to_edge_sq_distance<T>(da, db, dc));
        return distance;
    }

    template <>
    double DistanceBarrierRBProblem::distance<double>(
        const Eigen::VectorXd& sigma, const RB2Candidate& rbc)
    {
        typedef Eigen::Matrix<double, 3, 1> DTVector3d;
        typedef Eigen::Matrix<double, 2, 1> DTVector2d;

        DTVector3d sigma_E, sigma_V, position_E, position_V;

        sigma_V = sigma.segment(3 * rbc.vertex_body_id, 3);
        sigma_E = sigma.segment(3 * rbc.edge_body_id, 3);

        position_V = sigma_V.array()
            * m_assembler.m_dof_to_position.diagonal()
                  .segment(3 * rbc.vertex_body_id, 3)
                  .array();

        position_E = sigma_E.array()
            * m_assembler.m_dof_to_position.diagonal()
                  .segment(3 * rbc.edge_body_id, 3)
                  .array();

        const auto& rbs = m_assembler.m_rbs;
        DTVector2d da = rbs[size_t(rbc.edge_body_id)].world_vertex<double>(
            position_E, rbc.edge0_local_id);
        DTVector2d db = rbs[size_t(rbc.edge_body_id)].world_vertex<double>(
            position_E, rbc.edge1_local_id);
        DTVector2d dc = rbs[size_t(rbc.vertex_body_id)].world_vertex<double>(
            position_V, rbc.vertex_local_id);

        double distance = sqrt(point_to_edge_sq_distance<double>(da, db, dc));
        return distance;
    }

    void DistanceBarrierRBProblem::extract_local_system(
        const EdgeVertexCandidate& ev_candidate, RB2Candidate& rbc)
    {
        const long v_id = ev_candidate.vertex_index;
        const int e0_id = m_assembler.m_edges(ev_candidate.edge_index, 0);
        const int e1_id = m_assembler.m_edges(ev_candidate.edge_index, 1);
        int body_V_id, lv_id, body_E_id, le0_id, le1_id;

        m_assembler.global_to_local(v_id, body_V_id, lv_id);
        m_assembler.global_to_local(e0_id, body_E_id, le0_id);
        m_assembler.global_to_local(e1_id, body_E_id, le1_id);
        rbc.vertex_body_id = body_V_id;
        rbc.edge_body_id = body_E_id;
        rbc.vertex_local_id = lv_id;
        rbc.edge0_local_id = le0_id;
        rbc.edge1_local_id = le1_id;
    }

    Eigen::VectorXd DistanceBarrierRBProblem::compute_fd(
        const Eigen::VectorXd& sigma,
        const EdgeVertexCandidate& ev_candidate,
        const double h)
    {
        typedef AutodiffType<6> Diff;
        typedef Diff::DDouble1 T;
        Diff::activate();

        RB2Candidate rbc;
        extract_local_system(ev_candidate, rbc);
        T d = distance<T>(sigma, rbc);

        // distance finite diff
        auto f = [&](const Eigen::VectorXd& sigma_k) -> double {
            double dk = distance<double>(sigma_k, rbc);
            return dk;
        };

        Eigen::VectorXd approx_grad;
        finite_gradient(sigma, f, approx_grad, AccuracyOrder::SECOND, h);

        // chain rule
        double distance_grad = constraint_.distance_barrier_grad(d.getValue());
        approx_grad = approx_grad * distance_grad;
        return approx_grad;
    }

    bool DistanceBarrierRBProblem::compare_jac_g(const Eigen::VectorXd& sigma,
        const EdgeVertexCandidates& ev_candidates,
        const Eigen::MatrixXd& jac_g)
    {
        auto fx = eval_g(sigma);
        auto jac_full = eval_jac_g_full(sigma, ev_candidates);
        double norm = (jac_full - jac_g).norm();
        if (norm >= 1e-16) {
            spdlog::error("gradients dont match norm={:10e}", norm);
        }
        auto jac_approx = eval_jac_g_approx(sigma, ev_candidates);

        double norm_approx = (jac_full - jac_approx).norm();
        if (jac_full.size() != 0 && compare_jacobian(jac_approx, jac_full)) {
            spdlog::warn(
                "finite differences gradients dont match jac_diff_norm={:10e} gx.norm()={:.10e}",
                norm_approx, fx.norm());
            std::cout << sigma << std::endl;
        }
        return norm < 1e-16;
    }

    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g_approx(
        const Eigen::VectorXd& sigma, const EdgeVertexCandidates& ev_candidates)
    {

        auto func_g = [&](const Eigen::VectorXd& sigma_k) -> Eigen::VectorXd {
            Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma_k;
            Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

            EdgeVertexCandidates local_ev_candidates;
            auto check
                = constraint_.get_active_barrier_set(uk, local_ev_candidates);
            if (check != DistanceBarrierConstraint::NO_COLLISIONS) {
                throw std::logic_error("finite diff causes collision");
            }
            if (local_ev_candidates.size() != ev_candidates.size()) {
                throw std::logic_error("finite diff changes condidates set");
            }
            Eigen::VectorXd g_uk;
            constraint_.compute_candidates_constraints(uk, ev_candidates, g_uk);
            return g_uk;
        };

        Eigen::MatrixXd jac;
        bool success = false;
        double eps = 1e-6;
        while (!success && eps > 0.0) {
            try {
                ccd::finite_jacobian(
                    sigma, func_g, jac, AccuracyOrder::SECOND, eps);
                success = true;
                spdlog::warn(
                    "finite differences computed with eps={:.18e} epsilon={:.18e}",
                    eps, this->get_barrier_epsilon());
            } catch (std::logic_error e) {
                spdlog::warn(e.what());
                eps = eps / 2.0;
            }
        }
        if (!success) {
            spdlog::warn("failed to compute finite differences");
        }
        spdlog::trace("rb_problem chec_finite_diff BEGIN");
        compare_jacobian(jac_approx, jac_full, 1e-4,
            fmt::format(
                "check_finite_diff h=1E-7 barrier_eps={:3e} x=approx y=autodiff",
                constraint_.get_barrier_epsilon()));
        spdlog::trace("rb_problem chec_finite_diff END");
        return norm < 1e-16;
    }

} // namespace opt
} // namespace ccd
