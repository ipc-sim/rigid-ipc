#include "distance_barrier_rb_problem.hpp"

#include <iostream>
#include <utils/tensor.hpp>

#include <autodiff/finitediff.hpp>
#include <logger.hpp>
#include <profiler.hpp>

#include <constants.hpp>
#include <multiprecision.hpp>

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
        nlohmann::json inner_solver_settings
            = params[params["barrier_solver"]["inner_solver"]
                         .get<std::string>()];
        opt_solver_.inner_solver_settings(inner_solver_settings);
        RigidBodyProblem::settings(params["rigid_body_problem"]);
    }

    nlohmann::json DistanceBarrierRBProblem::state() const
    {
        nlohmann::json json = RigidBodyProblem::state();
        if (debug_min_distance_ < 0) {
            json["min_distance"] = nullptr;
        } else {
            json["min_distance"] = debug_min_distance_;
        }
        return json;
    }

    bool DistanceBarrierRBProblem::simulation_step(const double time_step)
    {
        bool has_collision = RigidBodyProblem::simulation_step(time_step);

        Eigen::VectorXd sigma
            = m_assembler.m_position_to_dof * m_assembler.rb_positions_t0();
        debug_min_distance_ = debug_min_distance(sigma);
        if (debug_min_distance_ >= 0) {
            spdlog::info(
                "candidate_step min_distance={:.8e}", debug_min_distance_);

            // our constraint is really d > min_d, we want to run
            // the optimization when we end the step with small
            // distances
            if (debug_min_distance_ <= constraint_.min_distance) {
                has_collision = true;
            }
        }

        return has_collision;
    }

    bool DistanceBarrierRBProblem::take_step(
        const Eigen::VectorXd& sigma, const double time_step)
    {
        debug_min_distance_ = debug_min_distance(sigma);
        if (debug_min_distance_ < 0) {
            spdlog::info("final_step min_distance=N/A");
        } else {
            spdlog::info("final_step min_distance={:.8e}", debug_min_distance_);
        }

        return RigidBodyProblem::take_step(sigma, time_step);
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

#if defined(DEBUG_LINESEARCH) || defined(DEBUG_COLLISIONS)
    Eigen::MatrixXd DistanceBarrierRBProblem::debug_vertices(
        const Eigen::VectorXd& sigma) const
    {
        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        return m_assembler.world_vertices(qk);
    }
#endif
    double DistanceBarrierRBProblem::debug_min_distance(
        const Eigen::VectorXd& sigma) const
    {
        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) - vertices_t0;

        Eigen::VectorXd d;
        constraint_.debug_compute_distances(uk, d);
        if (d.rows() > 0) {
            return d.minCoeff();
        }
        return -1;
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

    bool DistanceBarrierRBProblem::has_collisions(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j) const
    {
        Eigen::VectorXd qi = m_assembler.m_dof_to_position * sigma_i;
        Eigen::MatrixXd xi = m_assembler.world_vertices(qi);

        Eigen::VectorXd qj = m_assembler.m_dof_to_position * sigma_j;
        Eigen::MatrixXd xj = m_assembler.world_vertices(qj);
        return constraint_.has_active_collisions(xi, xj);
    }

    Eigen::Matrix<Multiprecision, Eigen::Dynamic, 1>
    DistanceBarrierRBProblem::eval_mp_g(const Eigen::VectorXd& /*sigma*/)
    {
        //        Eigen::VectorXd qk = m_assembler.m_dof_to_position * sigma;
        //        Eigen::MatrixXd uk = m_assembler.world_vertices(qk) -
        //        vertices_t0;

        Eigen::Matrix<Multiprecision, Eigen::Dynamic, 1> g_uk;
        //        EdgeVertexCandidates ev_candidates;
        //        auto check = constraint_.get_active_barrier_set(uk,
        //        ev_candidates);

        //        if (check == DistanceBarrierConstraint::HAS_COLLISION) {
        g_uk.resize(1);
        g_uk(0) = Multiprecision(std::numeric_limits<double>::infinity(), 256);

        //        } else {
        //            Eigen::Matrix<Multiprecision, Eigen::Dynamic,
        //            Eigen::Dynamic> uk_mp;

        //            uk_mp.resizeLike(uk);
        //            for (int i = 0; i < uk.size(); ++i) {
        //                uk_mp(i) = Multiprecision(uk(i), 256);
        //            }
        //            constraint_.compute_candidates_constraints<Multiprecision>(
        //                uk_mp, ev_candidates, g_uk);
        //        }
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
        constraint_.get_active_barrier_set(uk, ev_candidates);
        PROFILE_END(UPDATE)

        PROFILE_START(EVAL)
        Eigen::MatrixXd gx_jacobian = eval_jac_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL)

#ifdef WITH_DERIVATIVE_CHECK
        bool derivative_check = true;
        if (derivative_check) {
            compare_jac_g(sigma, ev_candidates, gx_jacobian);
        }
#endif
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
        constraint_.get_active_barrier_set(uk, ev_candidates);
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
        constraint_.get_active_barrier_set(uk, ev_candidates);
        PROFILE_END(UPDATE)

        constraint_.compute_candidates_constraints(uk, ev_candidates, gx);

        PROFILE_START(EVAL_GRAD)
        gx_jacobian = eval_jac_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL_GRAD)

        PROFILE_START(EVAL_HESS)
        gx_hessian = eval_hessian_g_core(sigma, ev_candidates);
        PROFILE_END(EVAL_HESS)

#ifdef WITH_DERIVATIVE_CHECK
        bool derivative_check = true;
        if (derivative_check) {

            assert(compare_jac_g(sigma, ev_candidates, gx_jacobian));
        }
#endif
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

        //        T distance = sqrt(point_to_edge_sq_distance<T>(da, db, dc));
        T distance = point_to_edge_distance<T>(da, db, dc);
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

        //        double distance = sqrt(point_to_edge_sq_distance<double>(da,
        //        db, dc));
        double distance = point_to_edge_distance<double>(da, db, dc);
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

    void DistanceBarrierRBProblem::compare_fd(const Eigen::VectorXd& sigma,
        const EdgeVertexCandidate& ev_candidate,
        const Eigen::VectorXd& grad)
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

        // distance finite diff
        Eigen::VectorXd approx_grad;
        Eigen::VectorXd exact_grad(sigma.rows());
        Eigen::VectorXd local_exact_grad = d.getGradient();
        exact_grad.setZero();
        exact_grad.segment(3 * rbc.vertex_body_id, 3)
            = local_exact_grad.segment(0, 3);
        exact_grad.segment(3 * rbc.edge_body_id, 3)
            = local_exact_grad.segment(3, 3);

        finite_gradient(sigma, f, approx_grad, AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);
        if (!compare_gradient(approx_grad, exact_grad,
                Constants::FINITE_DIFF_TEST,
                fmt::format(
                    "check_finite_diff DISTANCE barrier_eps={:3e} d={:3e}",
                    constraint_.get_barrier_epsilon(), d.getValue()))) {
        }

        // barrier finite diff - chain rule
        double distance_grad = constraint_.distance_barrier_grad(d.getValue());
        approx_grad = approx_grad * distance_grad;

        compare_gradient(approx_grad, grad, Constants::FINITE_DIFF_TEST,
            fmt::format("check_finite_diff BARRIER barrier_eps={:3e} d={:3e}",
                constraint_.get_barrier_epsilon(), d.getValue()));
    }

    bool DistanceBarrierRBProblem::compare_jac_g(const Eigen::VectorXd& sigma,
        const EdgeVertexCandidates& ev_candidates,
        const Eigen::MatrixXd& jac_g)
    {

        auto jac_full = eval_jac_g_full(sigma, ev_candidates);

        bool pass = compare_jacobian(
            jac_full, jac_g, /*test_eps=*/Constants::FULL_GRADIENT_TEST);
        if (!pass) {
            spdlog::error("autodiff_gradients_dont_match");
        }

        Eigen::MatrixXd jac_approx(jac_g.rows(), jac_g.cols());
        assert(jac_approx.rows() == int(ev_candidates.size()));
        for (size_t i = 0; i < ev_candidates.size(); ++i) {
            const auto& ev = ev_candidates[i];
            compare_fd(sigma, ev, jac_full.row(int(i)));
        }

        return pass;
    }

} // namespace opt
} // namespace ccd
