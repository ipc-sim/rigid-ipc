#include "distance_barrier_rb_problem.hpp"

// IPC Toolkit
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/friction/friction.hpp>

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
        , static_friction_speed_bound(1e-3)
        , friction_iterations(1)
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

        // Friction
        static_friction_speed_bound =
            params["friction_constraints"]["static_friction_speed_bound"];
        friction_iterations = params["friction_constraints"]["iterations"];

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
        opt_result = solve_constraints();
        _has_intersections = take_step(opt_result.x);
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
        int pos_ndof = physics::Pose<double>::dim_to_pos_ndof(dim());
        int rot_ndof = physics::Pose<double>::dim_to_rot_ndof(dim());

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
                    Diff::d2vars(0, pose.dof()));

                Diff::DDouble2 dExi =
                    compute_body_energy<Diff::DDouble2>(body, pose_diff);

                energies[i] = dExi.getValue();
                gradi = dExi.getGradient();
                Eigen::MatrixXX6d hessi = dExi.getHessian();

                // Project dense block to make assembled matrix PSD
                hessi.topLeftCorner(pos_ndof, pos_ndof) = Eigen::project_to_psd(
                    hessi.topLeftCorner(pos_ndof, pos_ndof));
                // NOTE: This is handled with Tikhonov regularization instead
                // hessi.bottomRightCorner(rot_ndof, rot_ndof) =
                //     Eigen::project_to_pd(
                //         hessi.bottomRightCorner(rot_ndof, rot_ndof));

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
                    Diff::d1vars(0, pose.dof()));

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
        // A large mass (e.g., from a large ground plane) can affect the
        // accuracy.
        double mass_Linf =
            m_assembler.m_rb_mass_matrix.diagonal().lpNorm<Eigen::Infinity>();
        double tol = std::max(1e-8 * mass_Linf, 1e-4);
        if (compute_grad) {
            Eigen::VectorXd grad_approx = eval_grad_energy_approx(*this, x);
            if (!fd::compare_gradient(grad, grad_approx, tol)) {
                spdlog::error("finite gradient check failed for E(x)");
            }
        }
        if (compute_hess) {
            Eigen::MatrixXd hess_approx = eval_hess_energy_approx(*this, x);
            if (!fd::compare_jacobian(hess, hess_approx, tol)) {
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

            // Eigen::Matrix3d Qdot_t0 = Q_t0 *
            // Eigen::Hat(body.velocity.rotation);
            Eigen::Matrix3d Qdot_t0 = body.Qdot;

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

        NAMED_PROFILE_POINT("compute_barrier__compute_constraints", COMPUTE)
        PROFILE_START(COMPUTE)
        double Bx = compute_barrier_term(
            x, candidates, grad, hess, compute_grad, compute_hess);
        PROFILE_END(COMPUTE)

#ifdef WITH_DERIVATIVE_CHECK
        if (compute_grad) {
            check_grad_barrier(x, candidates);
        }
        if (compute_hess) {
            check_hess_barrier(x, candidates);
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

        Eigen::VectorXd grad_f;
        Eigen::SparseMatrix<double> hess_f;
        double Fx = compute_friction_term(
            x, grad_f, hess_f, compute_grad, compute_hess);

        if (compute_grad) {
            grad += grad_f;
        }
        if (compute_hess) {
            hess += hess_f;
        }

        return Bx + Fx;
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
        NAMED_PROFILE_POINT("compute_barrier__compute_hess", COMPUTE_HESS)
        NAMED_PROFILE_POINT("compute_barrier__compute_grad", COMPUTE_GRAD)
        NAMED_PROFILE_POINT("compute_barrier__compute_val", COMPUTE_VAL)
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
            PROFILE_START(COMPUTE_GRAD)

            Diff::DDouble1 dBxi = distance_barrier<Diff::DDouble1>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();

            PROFILE_END(COMPUTE_GRAD)
        } else {
            PROFILE_START(COMPUTE_VAL)
            Bx += distance_barrier<double>(sigma, rbc);
            PROFILE_END(COMPUTE_VAL)
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

    // Compute the derivatives of the function that maps from rigid body DOF
    // to vertex positions.
    Eigen::MatrixXd DistanceBarrierRBProblem::rigid_dof_to_vertices(
        const Eigen::VectorXd& x,
        Eigen::MatrixXd& jac,
        std::vector<Eigen::MatrixXd>& hess,
        bool compute_jac,
        bool compute_hess)
    {
        // V: Rᵐ ↦ Rⁿ (vertices flattened rowwise)
        int m = x.size();
        int n = m_assembler.num_vertices() * dim();
        int rb_ndof = physics::Pose<double>::dim_to_ndof(dim());

        // Activate autodiff with the correct number of variables
        typedef AutodiffType<Eigen::Dynamic, /*maxN=*/6> Diff;
        Diff::activate(rb_ndof);

        if (compute_jac) {
            // ∇ V(x): Rᵐ ↦ Rⁿˣᵐ
            jac = Eigen::MatrixXd::Zero(n, m);
        }
        if (compute_hess) {
            // ∇²V(x): Rᵐ ↦ Rⁿˣᵐˣᵐ
            hess.reserve(n);
        }

        Eigen::MatrixXd V(m_assembler.num_vertices(), dim());
        for (int rb_i = 0; rb_i < m_assembler.num_bodies(); rb_i++) {
            const physics::RigidBody& rb = m_assembler.m_rbs[rb_i];
            // Index of ribid bodies first vertex in the global vertices
            long rb_v0_i = m_assembler.m_body_vertex_id[rb_i];

            if (compute_hess) {
                Diff::D2MatrixXd V_diff =
                    rb.world_vertices(physics::Pose<Diff::DDouble2>(
                        Diff::d2vars(0, x.segment(rb_i * rb_ndof, rb_ndof))));

                for (int i = 0; i < V_diff.rows(); i++) {
                    for (int j = 0; j < V_diff.cols(); j++) {
                        V(rb_v0_i + i, j) = V_diff(i, j).getValue();
                        jac.block(
                            /*i=*/(rb_v0_i + i) * dim() + j,
                            /*j=*/rb_i * rb_ndof, /*p=*/1, /*q=*/rb_ndof) =
                            V_diff(i, j).getGradient().transpose();
                        Eigen::MatrixXd full_hess = Eigen::MatrixXd::Zero(m, m);
                        full_hess.block(
                            rb_i * rb_ndof, rb_i * rb_ndof, rb_ndof, rb_ndof) =
                            V_diff(i, j).getHessian();
                        // hess[(rb_v0_i + i) * dim() + j] = full_hess;
                        hess.push_back(full_hess);
                    }
                }
            } else if (compute_jac) {
                Diff::D1MatrixXd V_diff =
                    rb.world_vertices(physics::Pose<Diff::DDouble1>(
                        Diff::d1vars(0, x.segment(rb_i * rb_ndof, rb_ndof))));

                for (int i = 0; i < V_diff.rows(); i++) {
                    for (int j = 0; j < V_diff.cols(); j++) {
                        V(rb_v0_i + i, j) = V_diff(i, j).getValue();
                        jac.block(
                            /*i=*/(rb_v0_i + i) * dim() + j,
                            /*j=*/rb_i * rb_ndof, /*p=*/1, /*q=*/rb_ndof) =
                            V_diff(i, j).getGradient().transpose();
                    }
                }

                return V;
            } else {
                V.block(rb_v0_i, 0, rb.vertices.rows(), rb.dim()) =
                    rb.world_vertices(physics::Pose<double>(
                        x.segment(rb_i * rb_ndof, rb_ndof)));
            }
        }

        return V;
    }

    double DistanceBarrierRBProblem::compute_friction_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        if (coefficient_friction == 0) {
            grad.setZero(x.size());
            hess = Eigen::SparseMatrix<double>(x.size(), x.size());
            return 0;
        }

        NAMED_PROFILE_POINT("compute_friction__update_constraints", UPDATE)
        PROFILE_START(UPDATE)
        physics::Poses<double> poses = this->dofs_to_poses(x);
        Candidates candidates;
        m_constraint.construct_active_barrier_set(
            m_assembler, poses, candidates);
        PROFILE_END(UPDATE)

        NAMED_PROFILE_POINT("compute_friction__compute_constraints", COMPUTE)
        PROFILE_START(COMPUTE)
        double friction_potential = compute_friction_term(
            x, candidates, grad, hess, compute_grad, compute_hess);
        PROFILE_END(COMPUTE)

#ifdef WITH_DERIVATIVE_CHECK
        if (compute_grad) {
            check_grad_friction(x, candidates, grad);
        }
        if (compute_hess) {
            check_hess_friction(x, candidates, hess);
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

        return friction_potential;
    }

    double DistanceBarrierRBProblem::compute_friction_term(
        const Eigen::VectorXd& x,
        const Candidates& distance_candidates,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);

        Eigen::MatrixXd jac_g;
        std::vector<Eigen::MatrixXd> hess_g;
        Eigen::MatrixXd V = rigid_dof_to_vertices(
            x, jac_g, hess_g, compute_grad || compute_hess, compute_hess);

        // Candidates to constraints
        ipc::Constraints constraints;
        for (const auto& ev_candidate : distance_candidates.ev_candidates) {
            constraints.ev_constraints.emplace_back(ev_candidate);
        }
        for (const auto& ee_candidate : distance_candidates.ee_candidates) {
            const auto& e00 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 0));
            const auto& e01 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 1));
            const auto& e10 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 0));
            const auto& e11 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 1));

            if (Eigen::cross(e01 - e00, e11 - e10).norm() > 1e-4) {
                double eps_x =
                    ipc::edge_edge_mollifier_threshold(e00, e01, e10, e11);
                constraints.ee_constraints.emplace_back(ee_candidate, eps_x);
            }
        }
        for (const auto& fv_candidate : distance_candidates.fv_candidates) {
            constraints.fv_constraints.emplace_back(fv_candidate);
        }

        // Contact constraints to friction constraints
        // TODO: Combine all these varaibles into the constraint set
        ipc::Constraints friction_constraint_set;
        std::vector<Eigen::VectorXd> closest_points;
        std::vector<Eigen::MatrixXd> tangent_bases;
        Eigen::VectorXd normal_force_magnitudes;
        ipc::compute_friction_bases(
            V0, m_assembler.m_edges, m_assembler.m_faces, constraints,
            barrier_activation_distance(), barrier_stiffness(),
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes);

        // Compute the surface friction potential
        double f = compute_friction_potential(
            V0, V, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes, static_friction_speed_bound * timestep(),
            coefficient_friction);

        Eigen::VectorXd grad_f;
        if (compute_grad || compute_hess) {
            // The gradient is also needed for the full hessian
            grad_f = compute_friction_potential_gradient(
                V0, V, m_assembler.m_edges, m_assembler.m_faces,
                friction_constraint_set, closest_points, tangent_bases,
                normal_force_magnitudes,
                static_friction_speed_bound * timestep(), coefficient_friction);
        }

        Eigen::SparseMatrix<double> hess_f;
        if (compute_hess) {
            hess_f = compute_friction_potential_hessian(
                V0, V, m_assembler.m_edges, m_assembler.m_faces,
                friction_constraint_set, closest_points, tangent_bases,
                normal_force_magnitudes,
                static_friction_speed_bound * timestep(), coefficient_friction);
        }

        // Apply the chain rule
        if (compute_grad) {
            // ∇ₓf(g(x)) = ∇ᵤf(u=g(x)) * ∇ₓg(x)
            grad = jac_g.transpose() * grad_f; // grad = [∇ₓf(g(x))]ᵀ
        }
        if (compute_hess) {
            Eigen::MatrixXd dense_hess = jac_g.transpose() * hess_f * jac_g;
            for (int i = 0; i < hess_g.size(); i++) {
                dense_hess += hess_g[i] * grad_f[i];
            }
            hess = dense_hess.sparseView();
        }

        return f;
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

    double DistanceBarrierRBProblem::compute_earliest_toi(
        const Eigen::VectorXd& sigma_i, const Eigen::VectorXd& sigma_j)
    {
        // m_use_barriers := solve_collisions
        // If we are not solve collisions then just compute if there was a
        // collision.
        if (!m_use_barriers) {
            this->has_collisions(sigma_i, sigma_j); // will set m_had_collisions
            return std::numeric_limits<double>::infinity();
        }

        physics::Poses<double> poses_i = this->dofs_to_poses(sigma_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(sigma_j);
        double earliest_toi =
            m_constraint.compute_earliest_toi(m_assembler, poses_i, poses_j);
        m_had_collisions |= earliest_toi <= 1;
        return earliest_toi;
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
            T c_div_eps_x = c / e_x;
            return (-c_div_eps_x + 2) * c_div_eps_x;
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
        const Eigen::VectorXd& sigma, const Candidates& candidates)
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
        const Eigen::VectorXd& sigma, const Candidates& candidates)
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

    void DistanceBarrierRBProblem::check_grad_friction(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::VectorXd& grad)
    {
        auto f = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_f;
            Eigen::SparseMatrix<double> hess_f;
            return compute_friction_term(
                x, candidates, grad_f, hess_f,
                /*compute_grad=*/false, /*compute_hess=*/false);
        };
        Eigen::VectorXd grad_approx;
        fd::finite_gradient(sigma, f, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx)) {
            spdlog::error("finite gradient check failed for friction");
        }

        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(sigma.size());
        Diff::D1MatrixXd V_diff = m_assembler.world_vertices(
            this->dofs_to_poses(Diff::d1vars(0, sigma)));

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);
        Eigen::MatrixXd V =
            m_assembler.world_vertices(this->dofs_to_poses(sigma));

        // Candidates to constraints
        ipc::Constraints constraints;
        for (const auto& ev_candidate : candidates.ev_candidates) {
            constraints.ev_constraints.emplace_back(ev_candidate);
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            const auto& e00 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 0));
            const auto& e01 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 1));
            const auto& e10 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 0));
            const auto& e11 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 1));

            if (Eigen::cross(e01 - e00, e11 - e10).norm() > 1e-4) {
                double eps_x =
                    ipc::edge_edge_mollifier_threshold(e00, e01, e10, e11);
                constraints.ee_constraints.emplace_back(ee_candidate, eps_x);
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            constraints.fv_constraints.emplace_back(fv_candidate);
        }

        // Contact constraints to friction constraints
        // TODO: Combine all these varaibles into the constraint set
        ipc::Constraints friction_constraint_set;
        std::vector<Eigen::VectorXd> closest_points;
        std::vector<Eigen::MatrixXd> tangent_bases;
        Eigen::VectorXd normal_force_magnitudes;
        ipc::compute_friction_bases(
            V0, m_assembler.m_edges, m_assembler.m_faces, constraints,
            barrier_activation_distance(), barrier_stiffness(),
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes);

        Diff::DDouble1 f_diff = compute_friction_potential(
            V0, V_diff, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes, static_friction_speed_bound * timestep(),
            coefficient_friction);

        if (std::isfinite(f_diff.getGradient().sum())
            && !fd::compare_gradient(grad, f_diff.getGradient())) {
            spdlog::error("autodiff gradient check failed for friction");
        }
    }

    void DistanceBarrierRBProblem::check_hess_friction(
        const Eigen::VectorXd& sigma,
        const Candidates& candidates,
        const Eigen::SparseMatrix<double>& hess)
    {
        // Finite differences breaks when the displacements are zero.
        if ((m_assembler.world_vertices(this->dofs_to_poses(sigma))
             - m_assembler.world_vertices(poses_t0))
                .template lpNorm<Eigen::Infinity>()
            == 0) {
            return;
        }

        Eigen::MatrixXd dense_hess(hess);

        auto f = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_f;
            Eigen::SparseMatrix<double> hess_f;
            return compute_friction_term(
                x, candidates, grad_f, hess_f,
                // /*compute_grad=*/true, /*compute_hess=*/false);
                /*compute_grad=*/false, /*compute_hess=*/false);
            // return grad_f;
        };
        Eigen::MatrixXd hess_approx;
        fd::finite_hessian(sigma, f, hess_approx);
        if (!fd::compare_hessian(hess, hess_approx)) {
            spdlog::error(
                "finite hessian check failed for friction "
                "(hess_L_inf_norm={:g} diff_L_inf_norm={:g})",
                dense_hess.lpNorm<Eigen::Infinity>(),
                (hess_approx - dense_hess).lpNorm<Eigen::Infinity>());
            std::cout << "hess_approx - dense_hess:\n"
                      << (1e-20 < (hess_approx - dense_hess).array().abs())
                             .select(hess_approx - dense_hess, 0.0f)
                      << std::endl;
        }
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(sigma.size());
        Diff::D2MatrixXd V_diff = m_assembler.world_vertices(
            this->dofs_to_poses(Diff::d2vars(0, sigma)));

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);
        Eigen::MatrixXd V =
            m_assembler.world_vertices(this->dofs_to_poses(sigma));

        // Candidates to constraints
        ipc::Constraints constraints;
        for (const auto& ev_candidate : candidates.ev_candidates) {
            constraints.ev_constraints.emplace_back(ev_candidate);
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            const auto& e00 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 0));
            const auto& e01 =
                V.row(m_assembler.m_edges(ee_candidate.edge0_index, 1));
            const auto& e10 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 0));
            const auto& e11 =
                V.row(m_assembler.m_edges(ee_candidate.edge1_index, 1));

            if (Eigen::cross(e01 - e00, e11 - e10).norm() > 1e-4) {
                double eps_x =
                    ipc::edge_edge_mollifier_threshold(e00, e01, e10, e11);
                constraints.ee_constraints.emplace_back(ee_candidate, eps_x);
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            constraints.fv_constraints.emplace_back(fv_candidate);
        }

        // Contact constraints to friction constraints
        // TODO: Combine all these varaibles into the constraint set
        ipc::Constraints friction_constraint_set;
        std::vector<Eigen::VectorXd> closest_points;
        std::vector<Eigen::MatrixXd> tangent_bases;
        Eigen::VectorXd normal_force_magnitudes;
        ipc::compute_friction_bases(
            V0, m_assembler.m_edges, m_assembler.m_faces, constraints,
            barrier_activation_distance(), barrier_stiffness(),
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes);

        Diff::DDouble2 f_diff = compute_friction_potential(
            V0, V_diff, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraint_set, closest_points, tangent_bases,
            normal_force_magnitudes, static_friction_speed_bound * timestep(),
            coefficient_friction);

        if (std::isfinite(f_diff.getHessian().sum())
            && !fd::compare_hessian(dense_hess, f_diff.getHessian())) {
            spdlog::error(
                "autodiff hessian check failed for friction "
                "(hess_L_inf_norm={:g} diff_L_inf_norm={:g})",
                dense_hess.lpNorm<Eigen::Infinity>(),
                (f_diff.getHessian() - dense_hess).lpNorm<Eigen::Infinity>());
            std::cout << "dense_hess - f_diff.getHessian():\n"
                      << (1e-20
                          < (dense_hess - f_diff.getHessian()).array().abs())
                             .select(dense_hess - f_diff.getHessian(), 0.0f)
                      << std::endl;
        } else if (!std::isfinite(f_diff.getHessian().sum())) {
            spdlog::warn("autodiff hessian failed for friction");
        }
    }
#endif

} // namespace opt
} // namespace ccd
