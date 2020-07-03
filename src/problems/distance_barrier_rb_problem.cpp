#include "distance_barrier_rb_problem.hpp"

#include <finitediff.hpp>

#include <constants.hpp>
#include <geometry/distance.hpp>
#include <solvers/solver_factory.hpp>
#include <utils/not_implemented_error.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

namespace opt {

    DistanceBarrierRBProblem::DistanceBarrierRBProblem()
        : m_barrier_stiffness(1)
        , min_distance(-1)
        , m_had_collisions(false)
    {
    }

    void DistanceBarrierRBProblem::settings(const nlohmann::json& params)
    {
        m_constraint.settings(params["distance_barrier_constraint"]);
        // Select the optimization solver
        std::string solver_name = params["solver"].get<std::string>();
        m_opt_solver = SolverFactory::factory().get_barrier_solver(solver_name);
        m_opt_solver->settings(params[solver_name]);
        m_opt_solver->set_problem(*this);
        if (m_opt_solver->has_inner_solver()) {
            m_opt_solver->inner_solver().settings(
                params[m_opt_solver->inner_solver().name()]);
        }
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

    ////////////////////////////////////////////////////////////
    // Rigid Body Problem

    void DistanceBarrierRBProblem::simulation_step(
        bool& had_collisions, bool& _has_intersections, bool solve_collisions)
    {
        // Advance the poses, but leave the current pose unchanged for now.
        tbb::parallel_for(size_t(0), m_assembler.num_bodies(), [&](size_t i) {
            m_assembler.m_rbs[i].pose_prev = m_assembler.m_rbs[i].pose;
        });

        update_dof(); // TODO: Update this code

        // Reset m_had_collision which will be filled in by has_collisions().
        m_had_collisions = false;

        // Disable barriers if solve_collision == false
        this->m_use_barriers = solve_collisions;

        update_constraint();
        opt::OptimizationResults result = solve_constraints();
        _has_intersections = take_step(result.x);
        had_collisions = m_had_collisions;
    }

    bool DistanceBarrierRBProblem::take_step(const Eigen::VectorXd& sigma)
    {
        min_distance = compute_min_distance(sigma);
        if (min_distance < 0) {
            spdlog::info("final_step min_distance=N/A");
        } else {
            spdlog::info("final_step min_distance={:.8e}", min_distance);
        }

        return RigidBodyProblem::take_step(sigma);
    }

    ////////////////////////////////////////////////////////////
    // Barrier Problem

    // Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    double DistanceBarrierRBProblem::compute_energy_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        typedef AutodiffType<Eigen::Dynamic, /*maxN=*/6> Diff;

        int ndof = physics::Pose<double>::dim_to_ndof(dim());

        Eigen::VectorXd energies(m_assembler.num_bodies());
        if (compute_grad) {
            grad.resize(x.size());
        }
        tbb::concurrent_vector<Eigen::Triplet<double>> hess_triplets;
        if (compute_hess) {
            // Hessian is a block diagonal with (ndof x ndof) blocks
            hess_triplets.reserve(m_assembler.num_bodies() * ndof * ndof);
        }

        const std::vector<physics::Pose<double>> poses = this->dofs_to_poses(x);
        assert(poses.size() == m_assembler.num_bodies());

        tbb::parallel_for(size_t(0), poses.size(), [&](size_t i) {
            // Activate autodiff with the correct number of variables
            Diff::activate(ndof);

            const physics::Pose<double>& pose = poses[i];
            const physics::RigidBody& body = m_assembler.m_rbs[i];

            Eigen::VectorX6d gradi;

            if (compute_hess) {
                // Initialize autodiff variables
                physics::Pose<Diff::DDouble2> pose_diff(
                    Diff::dTvars<Diff::DDouble2>(0, pose.dof()));

                Diff::DDouble2 dExi =
                    compute_body_energy<Diff::DDouble2>(body, pose_diff);

                energies[i] = dExi.getValue();
                gradi = dExi.getGradient();
                Eigen::MatrixXX6d hessi = dExi.getHessian();

                // Project dense block to make assembled matrix PSD
                // hessi = Eigen::project_to_psd(hessi);

                // Add the local hessian as triplets in the global hessian
                for (int r = 0; r < hessi.rows(); r++) {
                    for (int c = 0; c < hessi.cols(); c++) {
                        hess_triplets.emplace_back(
                            i * ndof + r, i * ndof + c, hessi(r, c));
                    }
                }
            } else if (compute_grad) {
                // Initialize autodiff variables
                physics::Pose<Diff::DDouble1> pose_diff(
                    Diff::dTvars<Diff::DDouble1>(0, pose.dof()));

                Diff::DDouble1 dExi =
                    compute_body_energy<Diff::DDouble1>(body, pose_diff);

                energies[i] = dExi.getValue();
                gradi = dExi.getGradient();
            } else {
                energies[i] = compute_body_energy<double>(body, pose);
            }

            if (compute_grad) {
                grad.segment(i * ndof, ndof) = gradi;
            }
        });

        if (compute_hess) {
            // ∇²E: Rⁿ ↦ Rⁿˣⁿ
            hess.resize(x.size(), x.size());
            hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
        }

#ifdef WITH_DERIVATIVE_CHECK
        if (compute_grad) {
            Eigen::VectorXd grad_approx = eval_grad_energy_approx(*this, x);
            if (!fd::compare_gradient(grad, grad_approx)) {
                spdlog::error("finite gradient check failed for E(x)");
            }
        }
        if (compute_hess) {
            Eigen::MatrixXd hess_approx = eval_hess_energy_approx(*this, x);
            if (!fd::compare_jacobian(hess, hess_approx)) {
                spdlog::error("finite hessian check failed for E(x)");
            }
        }
#endif

        return energies.sum();
    }

    // Compute the energy term for a single rigid body
    template <typename T>
    T DistanceBarrierRBProblem::compute_body_energy(
        const physics::RigidBody& body, const physics::Pose<T>& pose)
    {
        // NOTE: t0 suffix indicates the current value not the inital value
        double h = timestep();

        // Linear energy
        Eigen::VectorX3<T> q = pose.position;
        const Eigen::VectorX3d& q_t0 = body.pose.position;
        const Eigen::VectorX3d& qdot_t0 = body.velocity.position;
        Eigen::VectorX3d qdotdot_t0 = gravity + body.force.position / body.mass;

        // ½mqᵀq - mqᵀ(qᵗ + h(q̇ᵗ + h(g + f/m)))
        T energy = 0.5 * body.mass * q.dot(q)
            - body.mass
                * q.dot((q_t0 + h * (qdot_t0 + h * (qdotdot_t0))).cast<T>());

        // Rotational energy
        if (dim() == 3) {
            Eigen::Matrix3<T> Q = pose.construct_rotation_matrix();
            Eigen::Matrix3d Q_t0 = body.pose.construct_rotation_matrix();

            Eigen::Matrix3d Qdot_t0 = Q_t0 * Eigen::Hat(body.velocity.rotation);

            const Eigen::VectorX3d& I = body.moment_of_inertia;
            Eigen::DiagonalMatrix<T, 3> J(
                T(0.5 * (-I.x() + I.y() + I.z())),
                T(0.5 * (I.x() - I.y() + I.z())),
                T(0.5 * (I.x() + I.y() - I.z())));

            // ½tr(QJQᵀ) - tr(QJ(Qᵗ + hQ̇ᵗ)ᵀ)
            // TODO: Add torque
            energy += 0.5 * (Q * J * Q.transpose()).trace()
                - (Q * J * (Q_t0 + h * Qdot_t0).transpose().cast<T>()).trace();
        } else {
            assert(pose.rotation.size() == 1);
            // ½Iθ² - Iθ(θᵗ + hθ̇ᵗ)
            throw NotImplementedError(
                "DistanceBarrierRBProblem energy not implmented for 2D!");
        }

        return energy;
    }

    // Compute B(x) = ∑_{k ∈ C} b(d(x_k)) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
    double DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        int& num_constraints,
        bool compute_grad,
        bool compute_hess)
    {
        // Start by updating the constraint set
        NAMED_PROFILE_POINT("compute_barrier__update_constraints", UPDATE)
        PROFILE_START(UPDATE)
        physics::Poses<double> poses = this->dofs_to_poses(x);
        Candidates candidates;
        m_constraint.construct_active_barrier_set(
            m_assembler, poses, candidates);
        num_constraints = candidates.size();
        PROFILE_END(UPDATE)

        if (dim() == 2) {
            spdlog::debug(
                "problem={} num_edge_vertex_constraints={:d}", name(),
                candidates.ev_candidates.size());
        } else {
            spdlog::debug(
                "problem={} num_edge_edge_constraints={:d} "
                "num_face_vertex_constraints={:d}",
                name(), candidates.ee_candidates.size(),
                candidates.fv_candidates.size());
        }

        double Bx = compute_barrier_term(
            x, candidates, grad, hess, compute_grad, compute_hess);

#ifdef WITH_DERIVATIVE_CHECK
        if (compute_grad) {
            check_grad_barrier(x, candidates, grad);
        }
        if (compute_hess) {
            check_hess_barrier(x, candidates, grad);
        }
#endif

#ifndef NDEBUG
        // Make sure nothing went wrong
        if (compute_grad) {
            for (int i = 0; i < grad.size(); i++) {
                assert(std::isfinite(grad(i)));
            }
        }
        if (compute_hess) {
            typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;
            for (int k = 0; k < hess.outerSize(); ++k) {
                for (Iterator it(hess, k); it; ++it) {
                    assert(std::isfinite(it.value()));
                }
            }
        }
#endif
        return Bx;
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
        Eigen::VectorXd& grad,
        std::vector<Eigen::Triplet<double>>& hess_triplets,
        bool compute_grad,
        bool compute_hess)
    {
        // Activate autodiff with the correct number of variables
        typedef AutodiffType<Eigen::Dynamic, /*maxN=2*6=*/12> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyCandidate rbc;
        extract_local_system(candidate, rbc);
        std::array<long, 2> body_ids = rbc.body_ids();

        // Local gradient for a single constraint
        Eigen::VectorXd gradi;
        if (compute_hess) {
            NAMED_PROFILE_POINT("compute_barrier__compute_hess", COMPUTE_HESS)
            PROFILE_START(COMPUTE_HESS)

            Diff::DDouble2 dBxi = distance_barrier<Diff::DDouble2>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();
            Eigen::MatrixXd hessi = dBxi.getHessian();
            // Project dense block to make assembled matrix PSD
            hessi = Eigen::project_to_psd(hessi);
            // Add global triplets of the local values
            local_hessian_to_global_triplets(
                hessi, body_ids, ndof, hess_triplets);

            PROFILE_END(COMPUTE_HESS)
        } else if (compute_grad) {
            NAMED_PROFILE_POINT("compute_barrier__compute_grad", COMPUTE_GRAD)
            PROFILE_START(COMPUTE_GRAD)

            Diff::DDouble1 dBxi = distance_barrier<Diff::DDouble1>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();

            PROFILE_END(COMPUTE_GRAD)
        } else {
            Bx += distance_barrier<double>(sigma, rbc);
        }

        if (compute_grad) {
            assert(gradi.size() == 2 * ndof);
            grad.segment(ndof * body_ids[0], ndof) += gradi.head(ndof);
            grad.segment(ndof * body_ids[1], ndof) += gradi.tail(ndof);
        }
    }

    double DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& sigma,
        const Candidates& distance_candidates,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true)
    {
        // B: Rⁿ ↦ R
        double Bx = 0;
        // ∇B: Rⁿ ↦ Rⁿ
        grad = Eigen::VectorXd::Zero(num_vars_);
        // ∇²B: Rⁿ ↦ Rⁿˣⁿ
        std::vector<Eigen::Triplet<double>> hess_triplets;
        if (compute_hess) {
            int ndof = physics::Pose<double>::dim_to_ndof(dim());
            hess_triplets.reserve(distance_candidates.size() * 4 * ndof * ndof);
        }

        // Compute EV constraint hessian
        for (const auto& ev_candidate : distance_candidates.ev_candidates) {
            add_constraint_barrier<
                EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(
                sigma, ev_candidate, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        // Compute EE constraint hessian
        for (const auto& ee_candidate : distance_candidates.ee_candidates) {
            add_constraint_barrier<
                EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(
                sigma, ee_candidate, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        // Compute FV constraint hessian
        for (const auto& fv_candidate : distance_candidates.fv_candidates) {
            add_constraint_barrier<
                FaceVertexCandidate, RigidBodyFaceVertexCandidate>(
                sigma, fv_candidate, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        if (compute_hess) {
            // ∇²B: Rⁿ ↦ Rⁿˣⁿ
            hess.resize(num_vars_, num_vars_);
            hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
        }

        return Bx;
    }

    ///////////////////////////////////////////////////////////////////////////

    double DistanceBarrierRBProblem::compute_min_distance(
        const Eigen::VectorXd& sigma) const
    {
        physics::Poses<double> poses = this->dofs_to_poses(sigma);
        Eigen::VectorXd d;
        m_constraint.compute_distances(m_assembler, poses, d);
        if (d.rows() > 0) {
            return d.minCoeff();
        }
        return -1;
    }

    bool DistanceBarrierRBProblem::has_collisions(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j)
    {
        physics::Poses<double> poses_i = this->dofs_to_poses(sigma_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(sigma_j);
        bool collisions =
            m_constraint.has_active_collisions(m_assembler, poses_i, poses_j);
        m_had_collisions |= collisions;
        // m_use_barriers := solve_collisions
        return m_use_barriers ? collisions : false;
    }

    template <typename T, typename RigidBodyCandidate>
    T DistanceBarrierRBProblem::distance_barrier(
        const Eigen::VectorXd& sigma, const RigidBodyCandidate& rbc)
    {
        T d = distance<T>(sigma, rbc);
        T mollifier = constraint_mollifier<T>(sigma, rbc);
        T barrier = mollifier * m_constraint.distance_barrier<T>(d);
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
    T DistanceBarrierRBProblem::constraint_mollifier(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeEdgeCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_E0, sigma_E1;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.edge0_body_id, rbc.edge1_body_id, ndof, //
            sigma_E0, sigma_E1);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::Vector3<T> edge0_vertex0 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                sigma_E0, rbc.edge0_vertex0_local_id);
        Eigen::Vector3<T> edge0_vertex1 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                sigma_E0, rbc.edge0_vertex1_local_id);

        Eigen::Vector3<T> edge1_vertex0 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                sigma_E1, rbc.edge1_vertex0_local_id);
        Eigen::Vector3<T> edge1_vertex1 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                sigma_E1, rbc.edge1_vertex1_local_id);

        T c = (edge0_vertex1 - edge0_vertex0)
                  .cross(edge1_vertex1 - edge1_vertex0)
                  .squaredNorm();
        T e_x = 1e-3 * (edge0_vertex1 - edge0_vertex0).squaredNorm()
            * (edge1_vertex1 - edge1_vertex0).squaredNorm();
        if (c < e_x) {
            return -1 / (e_x * e_x) * c * c + 2 / e_x * c;
        } else {
            return T(1.0);
        }
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeVertexCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_V, sigma_E;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.vertex_body_id, rbc.edge_body_id, ndof, //
            sigma_V, sigma_E);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX<T> d_vertex = rbs[rbc.vertex_body_id].world_vertex<T>(
            sigma_V, rbc.vertex_local_id);
        Eigen::VectorX<T> d_edge_vertex0 =
            rbs[rbc.edge_body_id].world_vertex<T>(
                sigma_E, rbc.edge_vertex0_local_id);
        Eigen::VectorX<T> d_edge_vertex1 =
            rbs[rbc.edge_body_id].world_vertex<T>(
                sigma_E, rbc.edge_vertex1_local_id);

        // T distance = sqrt(point_to_edge_sq_distance<T>(da, db, dc));
        T distance = ccd::geometry::point_segment_distance<T>(
            d_vertex, d_edge_vertex0, d_edge_vertex1);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeEdgeCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_E0, sigma_E1;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.edge0_body_id, rbc.edge1_body_id, ndof, //
            sigma_E0, sigma_E1);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX<T> d_edge0_vertex0 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                sigma_E0, rbc.edge0_vertex0_local_id);
        Eigen::VectorX<T> d_edge0_vertex1 =
            rbs[rbc.edge0_body_id].world_vertex<T>(
                sigma_E0, rbc.edge0_vertex1_local_id);
        Eigen::VectorX<T> d_edge1_vertex0 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                sigma_E1, rbc.edge1_vertex0_local_id);
        Eigen::VectorX<T> d_edge1_vertex1 =
            rbs[rbc.edge1_body_id].world_vertex<T>(
                sigma_E1, rbc.edge1_vertex1_local_id);

        T distance = ccd::geometry::segment_segment_distance<T>(
            d_edge0_vertex0, d_edge0_vertex1, d_edge1_vertex0, d_edge1_vertex1);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyFaceVertexCandidate& rbc)
    {
        Eigen::VectorX6<T> sigma_V, sigma_F;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.vertex_body_id, rbc.face_body_id, ndof, //
            sigma_V, sigma_F);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX3<T> d_vertex = rbs[rbc.vertex_body_id].world_vertex<T>(
            sigma_V, rbc.vertex_local_id);
        Eigen::VectorX3<T> d_face_vertex0 =
            rbs[rbc.face_body_id].world_vertex<T>(
                sigma_F, rbc.face_vertex0_local_id);
        Eigen::VectorX3<T> d_face_vertex1 =
            rbs[rbc.face_body_id].world_vertex<T>(
                sigma_F, rbc.face_vertex1_local_id);
        Eigen::VectorX3<T> d_face_vertex2 =
            rbs[rbc.face_body_id].world_vertex<T>(
                sigma_F, rbc.face_vertex2_local_id);

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

    template <typename Candidate, typename RigidBodyCandidate>
    void DistanceBarrierRBProblem::check_distance_finite_gradient(
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

        std::array<long, 2> body_ids = rbc.body_ids();

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
                "check_finite_gradient DISTANCE barrier_eps={:3e} d={:3e}",
                m_constraint.barrier_activation_distance(), d.getValue()));
    }

    template <typename Candidate, typename RigidBodyCandidate>
    void DistanceBarrierRBProblem::check_distance_finite_hessian(
        const Eigen::VectorXd& sigma, const Candidate& candidate)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyCandidate rbc;
        extract_local_system(candidate, rbc);
        Diff::DDouble2 d = distance<Diff::DDouble2>(sigma, rbc);

        std::array<long, 2> body_ids = rbc.body_ids();

        // distance finite diff
        std::vector<Eigen::Triplet<double>> exact_hessian_triplets;
        local_hessian_to_global_triplets(
            d.getHessian(), body_ids, ndof, exact_hessian_triplets);
        Eigen::SparseMatrix<double> exact_hessian(sigma.rows(), sigma.rows());
        exact_hessian.setFromTriplets(
            exact_hessian_triplets.begin(), exact_hessian_triplets.end());

        // distance finite diff
        auto f = [&](const Eigen::VectorXd& sigma_k) {
            Eigen::VectorXd grad = Eigen::VectorXd::Zero(sigma_k.rows());
            Eigen::VectorXd local_grad =
                distance<Diff::DDouble1>(sigma_k, rbc).getGradient();
            grad.segment(ndof * body_ids[0], ndof) = local_grad.head(ndof);
            grad.segment(ndof * body_ids[1], ndof) = local_grad.tail(ndof);
            return grad;
        };
        Eigen::MatrixXd approx_hessian;
        finite_jacobian(
            sigma, f, approx_hessian, fd::AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);

        fd::compare_jacobian(
            approx_hessian, exact_hessian.toDense(),
            Constants::FINITE_DIFF_TEST,
            fmt::format(
                "check_finite_hessian DISTANCE barrier_eps={:3e} d={:3e}",
                m_constraint.barrier_activation_distance(), d.getValue()));
    }

    void DistanceBarrierRBProblem::check_grad_barrier(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::VectorXd& grad)
    {
        // Compute the finite difference of a single constraint
        for (const auto& ev : candidates.ev_candidates) {
            check_distance_finite_gradient<
                EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(sigma, ev);
        }
        for (const auto& ee : candidates.ee_candidates) {
            check_distance_finite_gradient<
                EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(sigma, ee);
        }
        for (const auto& fv : candidates.fv_candidates) {
            check_distance_finite_gradient<
                FaceVertexCandidate, RigidBodyFaceVertexCandidate>(sigma, fv);
        }
    }

    void DistanceBarrierRBProblem::check_hess_barrier(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::VectorXd& grad)
    {
        // Compute the finite difference of a single constraint
        for (const auto& ev : candidates.ev_candidates) {
            check_distance_finite_hessian<
                EdgeVertexCandidate, RigidBodyEdgeVertexCandidate>(sigma, ev);
        }
        for (const auto& ee : candidates.ee_candidates) {
            check_distance_finite_hessian<
                EdgeEdgeCandidate, RigidBodyEdgeEdgeCandidate>(sigma, ee);
        }
        for (const auto& fv : candidates.fv_candidates) {
            check_distance_finite_hessian<
                FaceVertexCandidate, RigidBodyFaceVertexCandidate>(sigma, fv);
        }
    }
#endif

} // namespace opt
} // namespace ccd
