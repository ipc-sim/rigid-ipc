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

    DistanceBarrierRBProblem::DistanceBarrierRBProblem()
        : min_distance(-1)
        , barrier_stiffness(1)
    {
    }

    void DistanceBarrierRBProblem::settings(const nlohmann::json& params)
    {
        constraint_.settings(params["distance_barrier_constraint"]);
        // TODO: Select the optimization solver
        opt_solver_.settings(params["barrier_solver"]);
        opt_solver_.set_problem(*this);
        nlohmann::json inner_solver_settings =
            params[params["barrier_solver"]["inner_solver"].get<std::string>()];
        opt_solver_.inner_solver().settings(inner_solver_settings);
        RigidBodyProblem::settings(params["rigid_body_problem"]);
    }

    nlohmann::json DistanceBarrierRBProblem::state() const
    {
        nlohmann::json json = RigidBodyProblem::state();
        if (min_distance < 0) {
            json["min_distance"] = nullptr;
        } else {
            json["min_distance"] = min_distance;
        }
        return json;
    }

    bool DistanceBarrierRBProblem::simulation_step(const double time_step)
    {
        // Take an unconstrained step
        bool has_collision = RigidBodyProblem::simulation_step(time_step);

        // Check if minimum distance is violated
        min_distance = compute_min_distance(
            this->poses_to_dofs(m_assembler.rb_poses_t0()));
        if (min_distance >= 0) {
            spdlog::info("candidate_step min_distance={:.8e}", min_distance);

            // our constraint is really d > min_d, we want to run the
            // optimization when we end the step with small distances
            if (min_distance <= constraint_.min_distance) {
                has_collision = true;
            }
        }

        return has_collision;
    }

    bool DistanceBarrierRBProblem::take_step(
        const Eigen::VectorXd& sigma, const double time_step)
    {
        min_distance = compute_min_distance(sigma);
        if (min_distance < 0) {
            spdlog::info("final_step min_distance=N/A");
        } else {
            spdlog::info("final_step min_distance={:.8e}", min_distance);
        }

        return RigidBodyProblem::take_step(sigma, time_step);
    }

    ////////////////////////////////////////////////////////////
    /// Barrier Problem

    // Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    void DistanceBarrierRBProblem::compute_energy_term(
        const Eigen::VectorXd& x,
        double& Ex,
        Eigen::VectorXd& grad_Ex,
        Eigen::SparseMatrix<double>& hess_Ex,
        bool compute_grad,
        bool compute_hess)
    {
        Eigen::VectorXd diff = x - this->poses_to_dofs(poses_t1);

        Eigen::SparseMatrix<double> M = m_assembler.m_dof_to_pose.transpose()
            * m_assembler.m_rb_mass_matrix * m_assembler.m_dof_to_pose;

        Ex = 0.5 * diff.transpose() * M * diff;

        if (compute_grad) {
            grad_Ex = M * diff;
#ifdef WITH_DERIVATIVE_CHECK
            Eigen::VectorXd grad_Ex_approx = eval_grad_energy_approx(*this, x);
            if (!fd::compare_gradient(grad_Ex, grad_Ex_approx)) {
                spdlog::error("finite gradient check failed for E(x)");
            }
#endif
        }

        if (compute_hess) {
            hess_Ex = M;
#ifdef WITH_DERIVATIVE_CHECK
            Eigen::MatrixXd hess_Ex_approx = eval_hess_energy_approx(*this, x);
            if (!fd::compare_jacobian(hess_Ex, hess_Ex_approx)) {
                spdlog::error("finite hessian check failed for E(x)");
            }
#endif
        }
    }

    // Compute B(x) = ∑_{k ∈ C} b(d(x_k)) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    int DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& x,
        double& Bx,
        Eigen::VectorXd& grad_Bx,
        Eigen::SparseMatrix<double>& hess_Bx,
        bool compute_grad,
        bool compute_hess)
    {
        // Start by updating the constraint set
        NAMED_PROFILE_POINT("compute_barrier__update_constraints", UPDATE)
        PROFILE_START(UPDATE)
        physics::Poses<double> poses = this->dofs_to_poses(x);
        Candidates candidates;
        constraint_.construct_active_barrier_set(
            m_assembler, poses, candidates);
        PROFILE_END(UPDATE)

        // if we want a derivative we use the autodiff version
        if (compute_grad || compute_hess) {
            compute_barrier_term(
                x, candidates, Bx, grad_Bx, hess_Bx, //
                compute_grad, compute_hess);
        } else {
            // TODO: This can be remove in favor of using the autodiff version
            // with double
            Eigen::VectorXd gx;
            constraint_.compute_candidates_constraints(
                m_assembler, poses, candidates, gx);
            Bx = gx.sum();
        }

#ifdef WITH_DERIVATIVE_CHECK
        if (compute_grad) {
            check_grad_barrier(x, candidates, grad_Bx);
        }
        // if (compute_hess) {
        //     check_hess_barrier(x, candidates, grad_Bx);
        // }
#endif

#ifndef NDEBUG
        // Make sure nothing went wrong
        if (compute_grad) {
            for (int i = 0; i < grad_Bx.size(); i++) {
                assert(std::isfinite(grad_Bx(i)));
            }
        }
        if (compute_hess) {
            typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;
            for (int k = 0; k < hess_Bx.outerSize(); ++k) {
                for (Iterator it(hess_Bx, k); it; ++it) {
                    assert(std::isfinite(it.value()));
                }
            }
        }
#endif

        // Return the number of constraints
        return candidates.size();
    }

    // Convert from a local hessian to the triplets in the global hessian
    void local_hessian_to_global_triplets(
        const Eigen::MatrixXd& local_hessian,
        const std::array<long, 2>& body_ids,
        int ndof,
        std::vector<Eigen::Triplet<double>>& triplets)
    {
        assert(local_hessian.rows() == 2 * ndof);
        assert(local_hessian.cols() == 2 * ndof);
        triplets.reserve(triplets.size() + local_hessian.size());
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

    // Compute the derivatives of a single constraint
    template <typename Candidate, typename RigidBodyCandidate>
    void DistanceBarrierRBProblem::add_constraint_barrier(
        const Eigen::VectorXd& sigma,
        const Candidate& candidate,
        double& Bx,
        Eigen::VectorXd& grad_Bx,
        std::vector<Eigen::Triplet<double>>& hess_Bx_triplets,
        bool compute_grad,
        bool compute_hess)
    {
        // Activate autodiff with the correct number of variables
        typedef AutodiffType<Eigen::Dynamic, /*maxN=2*6=*/12> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyCandidate rbc;
        extract_local_system(candidate, rbc);
        std::array<long, 2> body_ids = rbc.get_body_ids();

        // Local gradient for a single constraint
        Eigen::VectorXd grad_Bxi;
        if (compute_hess) {
            NAMED_PROFILE_POINT("compute_barrier__compute_hess", COMPUTE_HESS)
            PROFILE_START(COMPUTE_HESS)

            Diff::DDouble2 dBxi = distance_barrier<Diff::DDouble2>(sigma, rbc);
            Bx += dBxi.getValue();
            grad_Bxi = dBxi.getGradient();
            Eigen::MatrixXd hess_Bxi = dBxi.getHessian();

            // Add global triplets of the local values
            local_hessian_to_global_triplets(
                hess_Bxi, body_ids, ndof, hess_Bx_triplets);

            PROFILE_END(COMPUTE_HESS)
        } else if (compute_grad) {
            NAMED_PROFILE_POINT("compute_barrier__compute_grad", COMPUTE_GRAD)
            PROFILE_START(COMPUTE_GRAD)

            Diff::DDouble1 dBxi = distance_barrier<Diff::DDouble1>(sigma, rbc);
            Bx += dBxi.getValue();
            grad_Bxi = dBxi.getGradient();

            PROFILE_END(COMPUTE_GRAD)
        } else {
            Bx += distance_barrier<double>(sigma, rbc);
        }

        if (compute_grad) {
            assert(grad_Bxi.size() == 2 * ndof);
            grad_Bx.segment(ndof * body_ids[0], ndof) += grad_Bxi.head(ndof);
            grad_Bx.segment(ndof * body_ids[1], ndof) += grad_Bxi.tail(ndof);
        }
    }

    void DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& sigma,
        const Candidates& distance_candidates,
        double& Bx,
        Eigen::VectorXd& grad_Bx,
        Eigen::SparseMatrix<double>& hess_Bx,
        bool compute_grad = true,
        bool compute_hess = true)
    {
        // B: Rⁿ ↦ R
        Bx = 0;
        // ∇B: Rⁿ ↦ Rⁿ
        grad_Bx = Eigen::VectorXd::Zero(num_vars_);
        // ∇²B: Rⁿ ↦ Rⁿˣⁿ
        std::vector<Eigen::Triplet<double>> hess_Bx_triplets;

        // Compute EV constraint hessian
        for (const auto& ev_candidate : distance_candidates.ev_candidates) {
            add_constraint_barrier<
                EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(
                sigma, ev_candidate, Bx, grad_Bx, hess_Bx_triplets,
                compute_grad, compute_hess);
        }

        // Compute EE constraint hessian
        for (const auto& ee_candidate : distance_candidates.ee_candidates) {
            add_constraint_barrier<
                EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(
                sigma, ee_candidate, Bx, grad_Bx, hess_Bx_triplets,
                compute_grad, compute_hess);
        }

        // Compute FV constraint hessian
        for (const auto& fv_candidate : distance_candidates.fv_candidates) {
            add_constraint_barrier<
                FaceVertexCandidate, RigidBodyFaceVertexCandidate>(
                sigma, fv_candidate, Bx, grad_Bx, hess_Bx_triplets,
                compute_grad, compute_hess);
        }

        if (compute_hess) {
            // ∇²B: Rⁿ ↦ Rⁿˣⁿ
            hess_Bx.resize(num_vars_, num_vars_);
            hess_Bx.setFromTriplets(
                hess_Bx_triplets.begin(), hess_Bx_triplets.end());
        }
    }

    ///////////////////////////////////////////////////////////////////////////

    double DistanceBarrierRBProblem::compute_min_distance(
        const Eigen::VectorXd& sigma) const
    {
        physics::Poses<double> poses = this->dofs_to_poses(sigma);
        Eigen::VectorXd d;
        constraint_.compute_distances(m_assembler, poses, d);
        if (d.rows() > 0) {
            return d.minCoeff();
        }
        return -1;
    }

    bool DistanceBarrierRBProblem::has_collisions(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j) const
    {
        physics::Poses<double> poses_i = this->dofs_to_poses(sigma_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(sigma_j);
        return constraint_.has_active_collisions(m_assembler, poses_i, poses_j);
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

#ifdef WITH_DERIVATIVE_CHECK
    // The following functions are used exclusivly to check that the
    // gradient and hessian match a finite difference version.

    Eigen::VectorXd DistanceBarrierRBProblem::compute_full_grad_barrier(
        const Eigen::VectorXd& sigma, const Candidates& candidates)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(num_vars_);
        assert(sigma.size() == num_vars_);

        Diff::D1VectorXd d_sigma = Diff::d1vars(0, sigma);

        physics::Poses<Diff::DDouble1> d_poses = this->dofs_to_poses(d_sigma);
        Diff::D1VectorXd dBxi;
        constraint_.compute_candidates_constraints<Diff::DDouble1>(
            m_assembler, d_poses, candidates, dBxi);

        assert(candidates.size() == dBxi.rows());
        return dBxi.sum().getGradient();
    }

    template <typename Candidate, typename RigidBodyCandidate>
    void DistanceBarrierRBProblem::check_distance_finite_diff(
        const Eigen::VectorXd& sigma, const Candidate& candidate)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyCandidate rbc;
        extract_local_system(candidate, rbc);
        Diff::DDouble1 d = distance<Diff::DDouble1>(sigma, rbc);

        // distance finite diff
        auto f = [&](const Eigen::VectorXd& sigma_k) -> double {
            return distance<double>(sigma_k, rbc);
        };

        std::array<long, 2> body_ids = rbc.get_body_ids();

        // distance finite diff
        Eigen::VectorXd exact_grad = Eigen::VectorXd::Zero(sigma.rows());
        Eigen::VectorXd local_exact_grad = d.getGradient();
        exact_grad.segment(ndof * body_ids[0], ndof) =
            local_exact_grad.head(ndof);
        exact_grad.segment(ndof * body_ids[1], ndof) =
            local_exact_grad.tail(ndof);

        Eigen::VectorXd approx_grad;
        finite_gradient(
            sigma, f, approx_grad, fd::AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);
        fd::compare_gradient(
            approx_grad, exact_grad, Constants::FINITE_DIFF_TEST,
            fmt::format(
                "check_finite_diff DISTANCE barrier_eps={:3e} d={:3e}",
                constraint_.get_barrier_epsilon(), d.getValue()));
    }

    void DistanceBarrierRBProblem::check_grad_barrier(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::VectorXd& grad_Bx)
    {
        // Compute the gradient using full autodiff
        fd::compare_jacobian(
            compute_full_grad_barrier(sigma, candidates), grad_Bx,
            /*test_eps=*/Constants::FULL_GRADIENT_TEST,
            "autodiff gradients do not match ∇B(x)");

        // Compute the finite difference of a single constraint
        for (const auto& ev : candidates.ev_candidates) {
            check_distance_finite_diff<
                EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(sigma, ev);
        }
        for (const auto& ee : candidates.ee_candidates) {
            check_distance_finite_diff<
                EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(sigma, ee);
        }
        for (const auto& fv : candidates.fv_candidates) {
            check_distance_finite_diff<
                FaceVertexCandidate, RigidBodyFaceVertexCandidate>(sigma, fv);
        }
    }

#endif

} // namespace opt
} // namespace ccd
