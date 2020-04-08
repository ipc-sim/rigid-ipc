#include "distance_barrier_rb_problem.hpp"

#include <finitediff.hpp>

#include <constants.hpp>
#include <geometry/distance.hpp>
#include <multiprecision.hpp>
#include <utils/tensor.hpp>

#include <logger.hpp>
#include <profiler.hpp>

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
        nlohmann::json inner_solver_settings =
            params[params["barrier_solver"]["inner_solver"].get<std::string>()];
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
        // Take an unconstrained step
        bool has_collision = RigidBodyProblem::simulation_step(time_step);

        // Update the optimization dof (arclength for rotations)
        Eigen::VectorXd sigma = this->poses_to_dofs(m_assembler.rb_poses_t0());
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
    Eigen::MatrixXd
    DistanceBarrierRBProblem::debug_vertices(const Eigen::VectorXd& sigma) const
    {
        Eigen::VectorXd qk = m_assembler.m_dof_to_pose * sigma;
        return m_assembler.world_vertices(qk);
    }
#endif
    double DistanceBarrierRBProblem::debug_min_distance(
        const Eigen::VectorXd& sigma) const
    {
        physics::Poses<double> poses = this->dofs_to_poses(sigma);
        Eigen::VectorXd d;
        constraint_.debug_compute_distances(m_assembler, poses, d);
        if (d.rows() > 0) {
            return d.minCoeff();
        }
        return -1;
    }

    Eigen::VectorXd
    DistanceBarrierRBProblem::eval_g(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(
            m_assembler, this->dofs_to_poses(sigma), g_uk);
        return g_uk;
    }

    bool DistanceBarrierRBProblem::has_collisions(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j) const
    {
        physics::Poses<double> poses_i = this->dofs_to_poses(sigma_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(sigma_j);
        return constraint_.has_active_collisions(m_assembler, poses_i, poses_j);
    }

    Eigen::MatrixXd
    DistanceBarrierRBProblem::eval_jac_g(const Eigen::VectorXd& sigma)
    {
        NAMED_PROFILE_POINT("eval_jac_g__update_constraints", UPDATE)
        NAMED_PROFILE_POINT("eval_jac_g__eval_jac", EVAL)

        PROFILE_START(UPDATE)
        Candidates candidates;
        constraint_.construct_active_barrier_set(
            m_assembler, this->dofs_to_poses(sigma), candidates);
        PROFILE_END(UPDATE)

        PROFILE_START(EVAL)
        Eigen::MatrixXd gx_jacobian = eval_jac_g_core(sigma, candidates);
        PROFILE_END(EVAL)

#ifdef WITH_DERIVATIVE_CHECK
        compare_jac_g(sigma, candidates, gx_jacobian);
#endif
        return gx_jacobian;
    }

    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g_full(
        const Eigen::VectorXd& sigma, const Candidates& candidates)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(num_vars_);
        assert(sigma.size() == num_vars_);

        Diff::D1VectorXd d_sigma = Diff::d1vars(0, sigma);

        physics::Poses<Diff::DDouble1> d_poses_t0 =
            physics::cast<double, Diff::DDouble1>(poses_t0);
        physics::Poses<Diff::DDouble1> d_poses = this->dofs_to_poses(d_sigma);
        Diff::D1VectorXd d_g_uk;
        constraint_.compute_candidates_constraints<Diff::DDouble1>(
            m_assembler, d_poses, candidates, d_g_uk);

        assert(candidates.size() == d_g_uk.rows());
        return Diff::get_gradient(d_g_uk);
    }

    std::vector<Eigen::SparseMatrix<double>>
    DistanceBarrierRBProblem::eval_hessian_g(const Eigen::VectorXd& sigma)
    {
        NAMED_PROFILE_POINT("eval_hess_g__update_constraints", UPDATE)
        NAMED_PROFILE_POINT("eval_hess_g__eval", EVAL)

        PROFILE_START(UPDATE)
        Candidates candidates;
        constraint_.construct_active_barrier_set(
            m_assembler, this->dofs_to_poses(sigma), candidates);
        PROFILE_END(UPDATE)

        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        PROFILE_START(EVAL)
        gx_hessian = eval_hessian_g_core(sigma, candidates);
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
        physics::Poses<double> poses = this->dofs_to_poses(sigma);
        Candidates candidates;
        constraint_.construct_active_barrier_set(
            m_assembler, poses, candidates);
        PROFILE_END(UPDATE)

        constraint_.compute_candidates_constraints(
            m_assembler, poses, candidates, gx);

        PROFILE_START(EVAL_GRAD)
        gx_jacobian = eval_jac_g_core(sigma, candidates);
        PROFILE_END(EVAL_GRAD)

        PROFILE_START(EVAL_HESS)
        gx_hessian = eval_hessian_g_core(sigma, candidates);
        PROFILE_END(EVAL_HESS)

#ifdef WITH_DERIVATIVE_CHECK
        compare_jac_g(sigma, candidates, gx_jacobian);
#endif
    }

    Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g_core(
        const Eigen::VectorXd& sigma, const Candidates& distance_candidates)
    {
        // ∇²g: Rⁿ ↦ Rᵐˣⁿ
        Eigen::MatrixXd jac_g =
            Eigen::MatrixXd::Zero(distance_candidates.size(), num_vars_);

        typedef AutodiffType<Eigen::Dynamic, /*maxN=2*6=*/12> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        // TODO: Parallelize
        for (size_t i = 0; i < distance_candidates.ev_candidates.size(); ++i) {
            const auto& ev_candidate = distance_candidates.ev_candidates[i];

            RigidBodyEdgeVertexCandidate rbc;
            extract_local_system(ev_candidate, rbc);

            Eigen::VectorXd gradient =
                distance_barrier<Diff::DDouble1>(sigma, rbc).getGradient();

            size_t cstr_id = i;
            jac_g.block(int(cstr_id), ndof* rbc.vertex_body_id, 1, ndof) =
                gradient.head(ndof).transpose();
            jac_g.block(int(cstr_id), ndof* rbc.edge_body_id, 1, ndof) =
                gradient.tail(ndof).transpose();
        }

        for (size_t i = 0; i < distance_candidates.ee_candidates.size(); ++i) {
            const auto& ee_candidate = distance_candidates.ee_candidates[i];

            RigidBodyEdgeEdgeCandidate rbc;
            extract_local_system(ee_candidate, rbc);

            Eigen::VectorXd gradient =
                distance_barrier<Diff::DDouble1>(sigma, rbc).getGradient();

            size_t cstr_id = i + distance_candidates.ev_candidates.size();
            jac_g.block(int(cstr_id), ndof* rbc.edge0_body_id, 1, ndof) =
                gradient.head(ndof).transpose();
            jac_g.block(int(cstr_id), ndof* rbc.edge1_body_id, 1, ndof) =
                gradient.tail(ndof).transpose();
        }

        for (size_t i = 0; i < distance_candidates.fv_candidates.size(); ++i) {
            const auto& fv_candidate = distance_candidates.fv_candidates[i];

            RigidBodyFaceVertexCandidate rbc;
            extract_local_system(fv_candidate, rbc);

            Eigen::VectorXd gradient =
                distance_barrier<Diff::DDouble1>(sigma, rbc).getGradient();

            size_t cstr_id = i + distance_candidates.ev_candidates.size()
                + distance_candidates.ee_candidates.size();
            jac_g.block(int(cstr_id), ndof* rbc.vertex_body_id, 1, ndof) =
                gradient.head(ndof).transpose();
            jac_g.block(int(cstr_id), ndof* rbc.face_body_id, 1, ndof) =
                gradient.tail(ndof).transpose();
        }

#ifndef NDEBUG
        for (int i = 0; i < jac_g.size(); i++) {
            assert(std::isfinite(jac_g(i)));
        }
#endif

        return jac_g;
    }

    void local_hessian_to_global_triplets(
        const Eigen::MatrixXd& local_hessian,
        const std::array<long, 2>& body_ids,
        int ndof,
        std::vector<Eigen::Triplet<double>>& triplets)
    {
        triplets.clear();
        triplets.reserve(/*(2 * ndof) x (2 * ndof) = */ 4 * ndof * ndof);
        for (int b_i = 0; b_i < body_ids.size(); b_i++) {
            for (int b_j = 0; b_j < body_ids.size(); b_j++) {
                for (int dof_i = 0; dof_i < ndof; dof_i++) {
                    for (int dof_j = 0; dof_j < ndof; dof_j++) {
                        double v = local_hessian(
                            ndof * b_i + dof_i, ndof * b_j + dof_j);
                        int r = ndof * body_ids[b_i] + dof_i;
                        int c = ndof * body_ids[b_j] + dof_j;
                        triplets.emplace_back(r, c, v);
                    }
                }
            }
        }
    }

    std::vector<Eigen::SparseMatrix<double>>
    DistanceBarrierRBProblem::eval_hessian_g_core(
        const Eigen::VectorXd& sigma, const Candidates& distance_candidates)
    {
        // ∇²g: Rⁿ ↦ Rᵐˣᵐˣⁿ
        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        gx_hessian.resize(distance_candidates.size());

        typedef AutodiffType<Eigen::Dynamic, /*maxN=2*6=*/12> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        std::vector<Eigen::Triplet<double>> triplets;

        // Compute EV constraint hessian
        for (size_t i = 0; i < distance_candidates.ev_candidates.size(); ++i) {
            const auto& ev_candidate = distance_candidates.ev_candidates[i];

            RigidBodyEdgeVertexCandidate rbc;
            extract_local_system(ev_candidate, rbc);

            local_hessian_to_global_triplets(
                distance_barrier<Diff::DDouble2>(sigma, rbc).getHessian(),
                { { rbc.vertex_body_id, rbc.edge_body_id } }, ndof, triplets);

            Eigen::SparseMatrix<double> global_el_hessian(num_vars_, num_vars_);
            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());

            size_t cstr_id = i;
            gx_hessian[cstr_id] = global_el_hessian;
        }

        // Compute EE constraint hessian
        for (size_t i = 0; i < distance_candidates.ee_candidates.size(); ++i) {
            const auto& ee_candidate = distance_candidates.ee_candidates[i];

            RigidBodyEdgeEdgeCandidate rbc;
            extract_local_system(ee_candidate, rbc);

            local_hessian_to_global_triplets(
                distance_barrier<Diff::DDouble2>(sigma, rbc).getHessian(),
                { { rbc.edge0_body_id, rbc.edge1_body_id } }, ndof, triplets);

            Eigen::SparseMatrix<double> global_el_hessian(num_vars_, num_vars_);
            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());

            size_t cstr_id = i + distance_candidates.ev_candidates.size();
            gx_hessian[cstr_id] = global_el_hessian;
        }

        // Compute FV constraint hessian
        for (size_t i = 0; i < distance_candidates.fv_candidates.size(); ++i) {
            const auto& fv_candidate = distance_candidates.fv_candidates[i];

            RigidBodyFaceVertexCandidate rbc;
            extract_local_system(fv_candidate, rbc);

            local_hessian_to_global_triplets(
                distance_barrier<Diff::DDouble2>(sigma, rbc).getHessian(),
                { { rbc.vertex_body_id, rbc.face_body_id } }, ndof, triplets);

            Eigen::SparseMatrix<double> global_el_hessian(num_vars_, num_vars_);
            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());

            size_t cstr_id = i + distance_candidates.ev_candidates.size()
                + distance_candidates.ee_candidates.size();
            gx_hessian[cstr_id] = global_el_hessian;
        }

#ifndef NDEBUG
        // Make sure nothing went wrong
        for (const auto& Hi : gx_hessian) {
            for (int k = 0; k < Hi.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(Hi, k); it;
                     ++it) {
                    assert(std::isfinite(it.value()));
                }
            }
        }
#endif

        return gx_hessian;
    }

    template <typename T, typename RigidBodyCandidate>
    T DistanceBarrierRBProblem::distance_barrier(
        const Eigen::VectorXd& sigma, const RigidBodyCandidate& rbc)
    {
        T d = distance<T>(sigma, rbc);
        T barrier = constraint_.distance_barrier<T>(d);
        return barrier;
    }

    template <typename T>
    void init_body_sigmas(
        const Eigen::VectorXd& sigma,
        long body_id0,
        long body_id1,
        int ndof,
        Eigen::VectorX6<T>& sigma0,
        Eigen::VectorX6<T>& sigma1)
    {
        typedef AutodiffType<Eigen::Dynamic, /*maxN=2*6=*/12> Diff;
        Diff::activate(2 * ndof);
        sigma0 = Diff::dTvars<T>(0, sigma.segment(ndof * body_id0, ndof));
        sigma1 = Diff::dTvars<T>(ndof, sigma.segment(ndof * body_id1, ndof));
    }

    template <>
    void init_body_sigmas(
        const Eigen::VectorXd& sigma,
        long body_id0,
        long body_id1,
        int ndof,
        Eigen::VectorX6d& sigma0,
        Eigen::VectorX6d& sigma1)
    {
        sigma0 = sigma.segment(ndof * body_id0, ndof);
        sigma1 = sigma.segment(ndof * body_id1, ndof);
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeVertexCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_V, sigma_E, pose_V, pose_E;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.vertex_body_id, rbc.edge_body_id, ndof, //
            sigma_V, sigma_E);

        // Convert body dof (rotation in arclength) to a ndof pose vector
        pose_V = sigma_V.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.vertex_body_id, ndof)
                  .template cast<T>()
                  .array();
        pose_E = sigma_E.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.edge_body_id, ndof)
                  .template cast<T>()
                  .array();

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX<T> d_vertex = rbs[rbc.vertex_body_id].world_vertex<T>(
            pose_V, rbc.vertex_local_id);
        Eigen::VectorX<T> d_edge_vertex0 =
            rbs[rbc.edge_body_id].world_vertex<T>(
                pose_E, rbc.edge_vertex0_local_id);
        Eigen::VectorX<T> d_edge_vertex1 =
            rbs[rbc.edge_body_id].world_vertex<T>(
                pose_E, rbc.edge_vertex1_local_id);

        // T distance = sqrt(point_to_edge_sq_distance<T>(da, db, dc));
        T distance = ccd::geometry::point_segment_distance<T>(
            d_vertex, d_edge_vertex0, d_edge_vertex1);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeEdgeCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_E0, sigma_E1, pose_E0, pose_E1;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.edge0_body_id, rbc.edge1_body_id, ndof, //
            sigma_E0, sigma_E1);

        // Convert body dof (rotation in arclength) to a ndof pose vector
        pose_E0 = sigma_E0.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.edge0_body_id, ndof)
                  .template cast<T>()
                  .array();
        pose_E1 = sigma_E1.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.edge1_body_id, ndof)
                  .template cast<T>()
                  .array();

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX<T> d_edge0_vertex0 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                pose_E0, rbc.edge0_vertex0_local_id);
        Eigen::VectorX<T> d_edge0_vertex1 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                pose_E0, rbc.edge0_vertex1_local_id);
        Eigen::VectorX<T> d_edge1_vertex0 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                pose_E1, rbc.edge1_vertex0_local_id);
        Eigen::VectorX<T> d_edge1_vertex1 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                pose_E1, rbc.edge1_vertex1_local_id);

        T distance = ccd::geometry::segment_segment_distance<T>(
            d_edge0_vertex0, d_edge0_vertex1, d_edge1_vertex0, d_edge1_vertex1);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyFaceVertexCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_V, sigma_F, pose_V, pose_F;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.vertex_body_id, rbc.face_body_id, ndof, //
            sigma_V, sigma_F);

        // Convert body dof (rotation in arclength) to a ndof pose vector
        pose_V = sigma_V.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.vertex_body_id, ndof)
                  .template cast<T>()
                  .array();
        pose_F = sigma_F.array()
            * m_assembler.m_dof_to_pose.diagonal()
                  .segment(ndof * rbc.face_body_id, ndof)
                  .template cast<T>()
                  .array();

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX3<T> d_vertex = rbs[rbc.vertex_body_id].world_vertex<T>(
            pose_V, rbc.vertex_local_id);
        Eigen::VectorX3<T> d_face_vertex0 =
            rbs[rbc.face_body_id].world_vertex<T>(
                pose_F, rbc.face_vertex0_local_id);
        Eigen::VectorX3<T> d_face_vertex1 =
            rbs[rbc.face_body_id].world_vertex<T>(
                pose_F, rbc.face_vertex1_local_id);
        Eigen::VectorX3<T> d_face_vertex2 =
            rbs[rbc.face_body_id].world_vertex<T>(
                pose_F, rbc.face_vertex2_local_id);

        T distance = ccd::geometry::point_triangle_distance<T>(
            d_vertex, d_face_vertex0, d_face_vertex1, d_face_vertex2);
        return distance;
    }

    void DistanceBarrierRBProblem::extract_local_system(
        const EdgeVertexCandidate& ev_candidate,
        RigidBodyEdgeVertexCandidate& rbc)
    {
        const long v_id = ev_candidate.vertex_index;
        const long e_v0_id = m_assembler.m_edges(ev_candidate.edge_index, 0);
        const long e_v1_id = m_assembler.m_edges(ev_candidate.edge_index, 1);

        long body_V_id, lv_id, body_E_id, le_v0_id, le_v1_id;
        m_assembler.global_to_local_vertex(v_id, body_V_id, lv_id);
        m_assembler.global_to_local_vertex(e_v0_id, body_E_id, le_v0_id);
        m_assembler.global_to_local_vertex(e_v1_id, body_E_id, le_v1_id);

        rbc.vertex_body_id = body_V_id;
        rbc.edge_body_id = body_E_id;
        rbc.vertex_local_id = lv_id;
        rbc.edge_vertex0_local_id = le_v0_id;
        rbc.edge_vertex1_local_id = le_v1_id;
    }

    void DistanceBarrierRBProblem::extract_local_system(
        const EdgeEdgeCandidate& ee_candidate, RigidBodyEdgeEdgeCandidate& rbc)
    {
        const long e0_v0_id = m_assembler.m_edges(ee_candidate.edge0_index, 0);
        const long e0_v1_id = m_assembler.m_edges(ee_candidate.edge0_index, 1);
        const long e1_v0_id = m_assembler.m_edges(ee_candidate.edge1_index, 0);
        const long e1_v1_id = m_assembler.m_edges(ee_candidate.edge1_index, 1);

        long body_E0_id, le0_v0_id, le0_v1_id, body_E1_id, le1_v0_id, le1_v1_id;
        m_assembler.global_to_local_vertex(e0_v0_id, body_E0_id, le0_v0_id);
        m_assembler.global_to_local_vertex(e0_v1_id, body_E0_id, le0_v1_id);
        m_assembler.global_to_local_vertex(e1_v0_id, body_E1_id, le1_v0_id);
        m_assembler.global_to_local_vertex(e1_v1_id, body_E1_id, le1_v1_id);

        rbc.edge0_body_id = body_E0_id;
        rbc.edge1_body_id = body_E1_id;
        rbc.edge0_vertex0_local_id = le0_v0_id;
        rbc.edge0_vertex1_local_id = le0_v1_id;
        rbc.edge1_vertex0_local_id = le1_v0_id;
        rbc.edge1_vertex1_local_id = le1_v1_id;
    }

    void DistanceBarrierRBProblem::extract_local_system(
        const FaceVertexCandidate& fv_candidate,
        RigidBodyFaceVertexCandidate& rbc)
    {
        const long v_id = fv_candidate.vertex_index;
        const long f_v0_id = m_assembler.m_faces(fv_candidate.face_index, 0);
        const long f_v1_id = m_assembler.m_faces(fv_candidate.face_index, 1);
        const long f_v2_id = m_assembler.m_faces(fv_candidate.face_index, 2);

        long body_V_id, lv_id, body_F_id, lf_v0_id, lf_v1_id, lf_v2_id;
        m_assembler.global_to_local_vertex(v_id, body_V_id, lv_id);
        m_assembler.global_to_local_vertex(f_v0_id, body_F_id, lf_v0_id);
        m_assembler.global_to_local_vertex(f_v1_id, body_F_id, lf_v1_id);
        m_assembler.global_to_local_vertex(f_v2_id, body_F_id, lf_v2_id);

        rbc.vertex_body_id = body_V_id;
        rbc.face_body_id = body_F_id;
        rbc.vertex_local_id = lv_id;
        rbc.face_vertex0_local_id = lf_v0_id;
        rbc.face_vertex1_local_id = lf_v1_id;
        rbc.face_vertex2_local_id = lf_v2_id;
    }

    inline std::array<long, 2> get_body_ids(RigidBodyEdgeVertexCandidate rbc)
    {
        return { { rbc.vertex_body_id, rbc.edge_body_id } };
    }
    inline std::array<long, 2> get_body_ids(RigidBodyEdgeEdgeCandidate rbc)
    {
        return { { rbc.edge0_body_id, rbc.edge1_body_id } };
    }
    inline std::array<long, 2> get_body_ids(RigidBodyFaceVertexCandidate rbc)
    {
        return { { rbc.vertex_body_id, rbc.face_body_id } };
    }

    template <typename Candidate, typename RigidBodyCandidate>
    bool DistanceBarrierRBProblem::compare_fd(
        const Eigen::VectorXd& sigma,
        const Candidate& candidate,
        const Eigen::VectorXd& grad)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyCandidate rbc;
        extract_local_system(candidate, rbc);
        Diff::DDouble1 d = distance<Diff::DDouble1>(sigma, rbc);

        // distance finite diff
        auto f = [&](const Eigen::VectorXd& sigma_k) -> double {
            double dk = distance<double>(sigma_k, rbc);
            return dk;
        };

        std::array<long, 2> body_ids = get_body_ids(rbc);

        // distance finite diff
        Eigen::VectorXd approx_grad;
        Eigen::VectorXd exact_grad(sigma.rows());
        Eigen::VectorXd local_exact_grad = d.getGradient();
        exact_grad.setZero();
        exact_grad.segment(ndof * body_ids[0], ndof) =
            local_exact_grad.head(ndof);
        exact_grad.segment(ndof * body_ids[1], ndof) =
            local_exact_grad.tail(ndof);

        finite_gradient(
            sigma, f, approx_grad, fd::AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);
        return fd::compare_gradient(
            approx_grad, exact_grad, Constants::FINITE_DIFF_TEST,
            fmt::format(
                "check_finite_diff DISTANCE barrier_eps={:3e} d={:3e}",
                constraint_.get_barrier_epsilon(), d.getValue()));
    }

    bool DistanceBarrierRBProblem::compare_jac_g(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::MatrixXd& jac_g)
    {
        auto jac_full = eval_jac_g_full(sigma, candidates);

        bool pass = fd::compare_jacobian(
            jac_full, jac_g, /*test_eps=*/Constants::FULL_GRADIENT_TEST);
        if (!pass) {
            spdlog::error("autodiff gradients do not match");
        }

        for (size_t i = 0; i < candidates.ev_candidates.size(); ++i) {
            const auto& ev = candidates.ev_candidates[i];
            compare_fd<EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(
                sigma, ev, jac_full.row(int(i)));
            compare_fd<EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(
                sigma, ev, jac_g.row(int(i)));
        }
        for (size_t i = 0; i < candidates.ee_candidates.size(); ++i) {
            const auto& ee = candidates.ee_candidates[i];
            compare_fd<EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(
                sigma, ee, jac_full.row(int(i)));
            compare_fd<EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(
                sigma, ee, jac_g.row(int(i)));
        }
        for (size_t i = 0; i < candidates.fv_candidates.size(); ++i) {
            const auto& fv = candidates.fv_candidates[i];
            compare_fd<FaceVertexCandidate, RigidBodyFaceVertexCandidate>(
                sigma, fv, jac_full.row(int(i)));
            compare_fd<FaceVertexCandidate, RigidBodyFaceVertexCandidate>(
                sigma, fv, jac_g.row(int(i)));
        }

        assert(pass);
        return pass;
    }

} // namespace opt
} // namespace ccd
