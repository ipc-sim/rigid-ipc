#include "distance_barrier_rb_problem.hpp"

// IPC Toolkit
#include <ipc/ipc.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

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

    // Compute the objective function:
    // f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k)) + ∑_{k ∈ C} D(x_k)
    double DistanceBarrierRBProblem::compute_objective(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {

        // Compute rigid body energy term
        NAMED_PROFILE_POINT(
            "compute_objective__compute_energy_term", COMPUTE_ENERGY_TERM);
        PROFILE_START(COMPUTE_ENERGY_TERM);
        double Ex =
            compute_energy_term(x, grad, hess, compute_grad, compute_hess);
        PROFILE_END(COMPUTE_ENERGY_TERM);

        // The following is used to disable constraints if desired
        // (useful for testing).
        if (!m_use_barriers) {
            return Ex;
        }

        // Compute a common constraint set to use for contacts and friction
        // Start by updating the constraint set
        NAMED_PROFILE_POINT(
            "compute_objective__construct_constraint_set",
            CONSTRUCT_CONSTRAINTS);
        PROFILE_START(CONSTRUCT_CONSTRAINTS);
        ipc::Constraints constraints;
        m_constraint.construct_constraint_set(
            m_assembler, this->dofs_to_poses(x), constraints);
        PROFILE_END(CONSTRUCT_CONSTRAINTS);

        spdlog::debug(
            "problem={} num_vertex_vertex_constraint={:d} "
            "num_edge_vertex_constraints={:d} num_edge_edge_constraints={:d} "
            "num_face_vertex_constraints={:d}",
            name(), constraints.vv_constraints.size(),
            constraints.ev_constraints.size(),
            constraints.ee_constraints.size(),
            constraints.fv_constraints.size());

        Eigen::VectorXd grad_Bx;
        Eigen::SparseMatrix<double> hess_Bx;
        NAMED_PROFILE_POINT(
            "compute_objective__compute_barrier_term", COMPUTE_BARRIER_TERM);
        PROFILE_START(COMPUTE_BARRIER_TERM);
        double Bx = compute_barrier_term(
            x, constraints, grad_Bx, hess_Bx, compute_grad, compute_hess);
        PROFILE_END(COMPUTE_BARRIER_TERM);

        // D(x) is the friction potential (Equation 15 in the IPC paper)
        Eigen::VectorXd grad_Dx;
        Eigen::SparseMatrix<double> hess_Dx;
        NAMED_PROFILE_POINT(
            "compute_objective__compute_friction_term", COMPUTE_FRICTION_TERM);
        PROFILE_START(COMPUTE_FRICTION_TERM);
        double Dx = compute_friction_term(
            x, constraints, grad_Dx, hess_Dx, compute_grad, compute_hess);
        PROFILE_END(COMPUTE_FRICTION_TERM);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                // Energy gradient is checked in compute_energy_term()
                check_grad_barrier(x, constraints, grad_Bx);
                check_grad_friction(x, constraints, grad_Dx);
            }
            if (compute_hess) {
                // Energy hessian is checked in compute_energy_term()
                check_hess_barrier(x, constraints, hess_Bx);
                check_hess_friction(x, constraints, hess_Dx);
            }
            is_checking_derivative = false;
        }
#endif

        // Sum all the potentials
        double kappa = barrier_stiffness();
        if (compute_grad) {
            grad += kappa * grad_Bx + grad_Dx;
        }
        if (compute_hess) {
            hess += kappa * hess_Bx + hess_Dx;
        }
        return Ex + kappa * Bx + Dx;
    }

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
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            // A large mass (e.g., from a large ground plane) can affect the
            // accuracy.
            double mass_Linf = m_assembler.m_rb_mass_matrix.diagonal()
                                   .lpNorm<Eigen::Infinity>();
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
            is_checking_derivative = false;
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
        NAMED_PROFILE_POINT(
            "compute_barrier__construct_constraint_set", CONSTRUCT_CONSTRAINTS);
        PROFILE_START(CONSTRUCT_CONSTRAINTS);
        physics::Poses<double> poses = this->dofs_to_poses(x);
        ipc::Constraints constraints;
        m_constraint.construct_constraint_set(m_assembler, poses, constraints);
        num_constraints = constraints.num_constraints();
        PROFILE_END(CONSTRUCT_CONSTRAINTS);

        spdlog::debug(
            "problem={} num_vertex_vertex_constraint={:d} "
            "num_edge_vertex_constraints={:d} num_edge_edge_constraints={:d} "
            "num_face_vertex_constraints={:d}",
            name(), constraints.vv_constraints.size(),
            constraints.ev_constraints.size(),
            constraints.ee_constraints.size(),
            constraints.fv_constraints.size());

        NAMED_PROFILE_POINT("compute_barrier__compute_constraints", COMPUTE);
        PROFILE_START(COMPUTE);
        double Bx = compute_barrier_term(
            x, constraints, grad, hess, compute_grad, compute_hess);
        PROFILE_END(COMPUTE);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_grad_barrier(x, constraints, grad);
            }
            if (compute_hess) {
                check_hess_barrier(x, constraints, hess);
            }
            is_checking_derivative = false;
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
    template <typename Constraint, typename RigidBodyConstraint>
    void DistanceBarrierRBProblem::add_constraint_barrier(
        const Eigen::VectorXd& sigma,
        const Constraint& constraint,
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

        RigidBodyConstraint rbc(m_assembler, constraint);
        std::array<long, 2> body_ids = rbc.body_ids();

        // Local gradient for a single constraint
        Eigen::VectorXd gradi;
        if (compute_hess) {
            NAMED_PROFILE_POINT("compute_barrier__compute_hess", COMPUTE_HESS);
            PROFILE_START(COMPUTE_HESS);

            Diff::DDouble2 dBxi = distance_barrier<Diff::DDouble2>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();
            Eigen::MatrixXd hessi = dBxi.getHessian();
            // Project dense block to make assembled matrix PSD
            hessi = Eigen::project_to_psd(hessi);
            // Add global triplets of the local values
            local_hessian_to_global_triplets(
                hessi, body_ids, ndof, hess_triplets);

            PROFILE_END(COMPUTE_HESS);
        } else if (compute_grad) {
            NAMED_PROFILE_POINT("compute_barrier__compute_grad", COMPUTE_GRAD);
            PROFILE_START(COMPUTE_GRAD);

            Diff::DDouble1 dBxi = distance_barrier<Diff::DDouble1>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();

            PROFILE_END(COMPUTE_GRAD);
        } else {
            NAMED_PROFILE_POINT("compute_barrier__compute_val", COMPUTE_VAL);
            PROFILE_START(COMPUTE_VAL);
            Bx += distance_barrier<double>(sigma, rbc);
            PROFILE_END(COMPUTE_VAL);
        }

        if (compute_grad) {
            assert(gradi.size() == 2 * ndof);
            grad.segment(ndof * body_ids[0], ndof) += gradi.head(ndof);
            grad.segment(ndof * body_ids[1], ndof) += gradi.tail(ndof);
        }
    }

    double DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& x,
        const ipc::Constraints& constraints,
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
            hess_triplets.reserve(constraints.size() * 4 * ndof * ndof);
        }

        // Compute VV constraint
        for (const auto& vv_constraint : constraints.vv_constraints) {
            add_constraint_barrier<
                ipc::VertexVertexConstraint, RigidBodyVertexVertexConstraint>(
                x, vv_constraint, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        // Compute EV constraint
        for (const auto& ev_constraint : constraints.ev_constraints) {
            add_constraint_barrier<
                ipc::EdgeVertexConstraint, RigidBodyEdgeVertexConstraint>(
                x, ev_constraint, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        // Compute EE constraint
        for (const auto& ee_constraint : constraints.ee_constraints) {
            add_constraint_barrier<
                ipc::EdgeEdgeConstraint, RigidBodyEdgeEdgeConstraint>(
                x, ee_constraint, Bx, grad, hess_triplets, compute_grad,
                compute_hess);
        }

        // Compute FV constraint
        for (const auto& fv_constraint : constraints.fv_constraints) {
            add_constraint_barrier<
                ipc::FaceVertexConstraint, RigidBodyFaceVertexConstraint>(
                x, fv_constraint, Bx, grad, hess_triplets, compute_grad,
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
            } else {
                V.block(rb_v0_i, 0, rb.vertices.rows(), rb.dim()) =
                    rb.world_vertices(physics::Pose<double>(
                        x.segment(rb_i * rb_ndof, rb_ndof)));
            }
        }

        assert(
            (V - m_assembler.world_vertices(this->dofs_to_poses(x))).norm()
            < 1e-12);

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
        ipc::Constraints constraints;
        m_constraint.construct_constraint_set(m_assembler, poses, constraints);
        PROFILE_END(UPDATE)

        NAMED_PROFILE_POINT("compute_friction__compute_constraints", COMPUTE)
        PROFILE_START(COMPUTE)
        double friction_potential = compute_friction_term(
            x, constraints, grad, hess, compute_grad, compute_hess);
        PROFILE_END(COMPUTE)

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_grad_friction(x, constraints, grad);
            }
            if (compute_hess) {
                check_hess_friction(x, constraints, hess);
            }
            is_checking_derivative = false;
        }
#endif

        return friction_potential;
    }

    double DistanceBarrierRBProblem::compute_friction_term(
        const Eigen::VectorXd& x,
        const ipc::Constraints& contact_constraints,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        if (contact_constraints.size() == 0 && coefficient_friction == 0) {
            grad.setZero(x.size());
            hess = Eigen::SparseMatrix<double>(x.size(), x.size());
            return 0;
        }

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);

        // Contact constraints to friction constraints
        ipc::FrictionConstraints friction_constraints;
        ipc::construct_friction_constraint_set(
            V0, m_assembler.m_edges, m_assembler.m_faces, contact_constraints,
            barrier_activation_distance(), barrier_stiffness(),
            coefficient_friction, friction_constraints);

        // Compute g(x)
        Eigen::MatrixXd jac_g;
        std::vector<Eigen::MatrixXd> hess_g;
        Eigen::MatrixXd V = rigid_dof_to_vertices(
            x, jac_g, hess_g, compute_grad || compute_hess, compute_hess);

        // Compute the surface friction potential
        double f = compute_friction_potential(
            V0, V, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraints, static_friction_speed_bound * timestep());

        Eigen::VectorXd grad_f;
        if (compute_grad || compute_hess) {
            // The gradient is also needed for the full hessian
            grad_f = compute_friction_potential_gradient(
                V0, V, m_assembler.m_edges, m_assembler.m_faces,
                friction_constraints, static_friction_speed_bound * timestep());
        }

        Eigen::SparseMatrix<double> hess_f;
        if (compute_hess) {
            hess_f = compute_friction_potential_hessian(
                V0, V, m_assembler.m_edges, m_assembler.m_faces,
                friction_constraints, static_friction_speed_bound * timestep(),
                /*project_to_psd=*/false);
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
            hess = Eigen::project_to_psd(dense_hess).sparseView();
        }

        return f;
    }

    ///////////////////////////////////////////////////////////////////////////

    double DistanceBarrierRBProblem::compute_min_distance(
        const Eigen::VectorXd& sigma) const
    {
        physics::Poses<double> poses = this->dofs_to_poses(sigma);
        double min_distance =
            m_constraint.compute_minimum_distance(m_assembler, poses);
        return std::isfinite(min_distance) ? min_distance : -1;
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

    template <typename T, typename RigidBodyConstraint>
    T DistanceBarrierRBProblem::distance_barrier(
        const Eigen::VectorXd& sigma, const RigidBodyConstraint& rbc)
    {
        return rbc.multiplicity * constraint_mollifier<T>(sigma, rbc)
            * m_constraint.distance_barrier<T>(distance<T>(sigma, rbc));
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
        const Eigen::VectorXd& sigma, const RigidBodyEdgeEdgeConstraint& rbc)
    {
        Eigen::VectorX6<T> sigma_E0, sigma_E1;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.edge0_body_id, rbc.edge1_body_id, ndof, //
            sigma_E0, sigma_E1);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::Vector3<T> ea0 = rbs[rbc.edge0_body_id].world_vertex<T>(
            sigma_E0, rbc.edge0_vertex0_local_id);
        Eigen::Vector3<T> ea1 = rbs[rbc.edge0_body_id].world_vertex<T>(
            sigma_E0, rbc.edge0_vertex1_local_id);

        Eigen::Vector3<T> eb0 = rbs[rbc.edge1_body_id].world_vertex<T>(
            sigma_E1, rbc.edge1_vertex0_local_id);
        Eigen::Vector3<T> eb1 = rbs[rbc.edge1_body_id].world_vertex<T>(
            sigma_E1, rbc.edge1_vertex1_local_id);

        return ipc::edge_edge_mollifier(ea0, ea1, eb0, eb1, rbc.eps_x);
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma,
        const RigidBodyVertexVertexConstraint& rbc)
    {
        Eigen::VectorX6<T> sigma_V0, sigma_V1;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        init_body_sigmas(
            sigma, rbc.vertex0_body_id, rbc.vertex1_body_id, ndof, //
            sigma_V0, sigma_V1);

        const auto& rbs = m_assembler.m_rbs;
        Eigen::VectorX<T> d_vertex0 = rbs[rbc.vertex0_body_id].world_vertex<T>(
            sigma_V0, rbc.vertex0_local_id);
        Eigen::VectorX<T> d_vertex1 = rbs[rbc.vertex1_body_id].world_vertex<T>(
            sigma_V1, rbc.vertex1_local_id);

        T distance = ipc::point_point_distance(d_vertex0, d_vertex1);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeVertexConstraint& rbc)
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
        T distance = ipc::point_edge_distance(
            d_vertex, d_edge_vertex0, d_edge_vertex1,
            ipc::PointEdgeDistanceType::P_E);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyEdgeEdgeConstraint& rbc)
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

        T distance = ipc::edge_edge_distance(
            d_edge0_vertex0, d_edge0_vertex1, d_edge1_vertex0, d_edge1_vertex1,
            ipc::EdgeEdgeDistanceType::EA_EB);
        return distance;
    }

    template <typename T>
    T DistanceBarrierRBProblem::distance(
        const Eigen::VectorXd& sigma, const RigidBodyFaceVertexConstraint& rbc)
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

        T distance = ipc::point_triangle_distance(
            d_vertex, d_face_vertex0, d_face_vertex1, d_face_vertex2,
            ipc::PointTriangleDistanceType::P_T);
        return distance;
    }

#ifdef WITH_DERIVATIVE_CHECK
    // The following functions are used exclusivly to check that the
    // gradient and hessian match a finite difference version.

    template <typename Constraint, typename RigidBodyConstraint>
    void DistanceBarrierRBProblem::check_distance_finite_gradient(
        const Eigen::VectorXd& sigma, const Constraint& constraint)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyConstraint rbc(m_assembler, constraint);
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
        fd::finite_gradient(
            sigma, f, approx_grad, fd::AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);
        if (!fd::compare_gradient(
                approx_grad, exact_grad, Constants::FINITE_DIFF_TEST)) {
            spdlog::error(
                "check_finite_gradient DISTANCE barrier_eps={:3e} d={:3e}",
                m_constraint.barrier_activation_distance(), d.getValue());
        }
    }

    template <typename Constraint, typename RigidBodyConstraint>
    void DistanceBarrierRBProblem::check_distance_finite_hessian(
        const Eigen::VectorXd& sigma, const Constraint& constraint)
    {
        typedef AutodiffType<Eigen::Dynamic> Diff;
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        Diff::activate(2 * ndof);

        RigidBodyConstraint rbc(m_assembler, constraint);
        Diff::DDouble2 d = distance<Diff::DDouble2>(sigma, rbc);

        std::array<long, 2> body_ids = rbc.body_ids();

        // distance autodiff
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
        fd::finite_jacobian(
            sigma, f, approx_hessian, fd::AccuracyOrder::SECOND,
            Constants::FINITE_DIFF_H);

        if (!fd::compare_jacobian(
                approx_hessian, exact_hessian.toDense(),
                Constants::FINITE_DIFF_TEST)) {
            spdlog::error(
                "check_finite_hessian DISTANCE barrier_eps={:3e} d={:3e}",
                m_constraint.barrier_activation_distance(), d.getValue());
        }
    }

    void DistanceBarrierRBProblem::check_grad_barrier(
        const Eigen::VectorXd& sigma,
        const ipc::Constraints& constraints,
        const Eigen::VectorXd& grad)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        for (int i = 0; i < grad.size(); i++) {
            if (!std::isfinite(grad(i))) {
                spdlog::error("barrier gradient is not finite");
            }
        }

        // Compute the finite difference of a single constraint
        for (const auto& vv : constraints.vv_constraints) {
            check_distance_finite_gradient<
                ipc::VertexVertexConstraint, RigidBodyVertexVertexConstraint>(
                sigma, vv);
        }
        for (const auto& ev : constraints.ev_constraints) {
            check_distance_finite_gradient<
                ipc::EdgeVertexConstraint, RigidBodyEdgeVertexConstraint>(
                sigma, ev);
        }
        for (const auto& ee : constraints.ee_constraints) {
            check_distance_finite_gradient<
                ipc::EdgeEdgeConstraint, RigidBodyEdgeEdgeConstraint>(
                sigma, ee);
        }
        for (const auto& fv : constraints.fv_constraints) {
            check_distance_finite_gradient<
                ipc::FaceVertexConstraint, RigidBodyFaceVertexConstraint>(
                sigma, fv);
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        auto b = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_b;
            Eigen::SparseMatrix<double> hess_b;
            return compute_barrier_term(
                x, constraints, grad_b, hess_b,
                /*compute_grad=*/false, /*compute_hess=*/false);
        };
        Eigen::VectorXd grad_approx;
        fd::finite_gradient(sigma, b, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx, 1e-3)) {
            spdlog::error("finite gradient check failed for barrier");
        }
    }

    void DistanceBarrierRBProblem::check_hess_barrier(
        const Eigen::VectorXd& sigma,
        const ipc::Constraints& constraints,
        const Eigen::SparseMatrix<double>& hess)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;
        for (int k = 0; k < hess.outerSize(); ++k) {
            for (Iterator it(hess, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    spdlog::error("barrier hessian is not finite");
                    return;
                }
            }
        }

        // Compute the finite difference of a single constraint
        for (const auto& vv : constraints.vv_constraints) {
            check_distance_finite_hessian<
                ipc::VertexVertexConstraint, RigidBodyVertexVertexConstraint>(
                sigma, vv);
        }
        for (const auto& ev : constraints.ev_constraints) {
            check_distance_finite_hessian<
                ipc::EdgeVertexConstraint, RigidBodyEdgeVertexConstraint>(
                sigma, ev);
        }
        for (const auto& ee : constraints.ee_constraints) {
            check_distance_finite_hessian<
                ipc::EdgeEdgeConstraint, RigidBodyEdgeEdgeConstraint>(
                sigma, ee);
        }
        for (const auto& fv : constraints.fv_constraints) {
            check_distance_finite_hessian<
                ipc::FaceVertexConstraint, RigidBodyFaceVertexConstraint>(
                sigma, fv);
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        // WARNING: The following check does not work well because the
        // different projections to PSD can affect results.
        // Eigen::MatrixXd dense_hess(hess);
        // auto f = [&](const Eigen::VectorXd& x) {
        //     Eigen::VectorXd grad_b;
        //     Eigen::SparseMatrix<double> hess_b;
        //     compute_barrier_term(
        //         x, constraints, grad_b, hess_b,
        //         /*compute_grad=*/true, /*compute_hess=*/false);
        //     return grad_b;
        // };
        // Eigen::MatrixXd hess_approx;
        // fd::finite_jacobian(sigma, f, hess_approx);
        // hess_approx = Eigen::project_to_psd(hess_approx);
        // if (!fd::compare_jacobian(
        //         hess, hess_approx, Constants::FINITE_DIFF_TEST)) {
        //     spdlog::error("finite hessian check failed for barrier");
        // }
    }

    void DistanceBarrierRBProblem::check_grad_friction(
        const Eigen::VectorXd& sigma,
        const ipc::Constraints& constraints,
        const Eigen::VectorXd& grad)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        for (int i = 0; i < grad.size(); i++) {
            if (!std::isfinite(grad(i))) {
                spdlog::error("friction gradient is not finite");
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        auto f = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_f;
            Eigen::SparseMatrix<double> hess_f;
            return compute_friction_term(
                x, constraints, grad_f, hess_f,
                /*compute_grad=*/false, /*compute_hess=*/false);
        };
        Eigen::VectorXd grad_approx;
        fd::finite_gradient(sigma, f, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx)) {
            spdlog::error("finite gradient check failed for friction");
        }

        ///////////////////////////////////////////////////////////////////////
        // Auto. diff. check
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(sigma.size());
        Diff::D1MatrixXd V_diff = m_assembler.world_vertices(
            this->dofs_to_poses(Diff::d1vars(0, sigma)));

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);

        // Contact constraints to friction constraints
        ipc::FrictionConstraints friction_constraint_set;
        ipc::construct_friction_constraint_set(
            V0, m_assembler.m_edges, m_assembler.m_faces, constraints,
            barrier_activation_distance(), barrier_stiffness(),
            coefficient_friction, friction_constraint_set);

        Diff::DDouble1 f_diff = compute_friction_potential(
            V0, V_diff, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraint_set, static_friction_speed_bound * timestep());

        if (std::isfinite(f_diff.getGradient().sum())
            && !fd::compare_gradient(grad, f_diff.getGradient())) {
            spdlog::error("autodiff gradient check failed for friction");
        }
    }

    void DistanceBarrierRBProblem::check_hess_friction(
        const Eigen::VectorXd& sigma,
        const ipc::Constraints& constraints,
        const Eigen::SparseMatrix<double>& hess)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;
        for (int k = 0; k < hess.outerSize(); ++k) {
            for (Iterator it(hess, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    spdlog::error("barrier hessian is not finite");
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        // Finite differences breaks when the displacements are zero.
        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);
        Eigen::MatrixXd V1 =
            m_assembler.world_vertices(this->dofs_to_poses(sigma));
        if ((V1 - V0).lpNorm<Eigen::Infinity>() == 0) {
            return;
        }

        Eigen::MatrixXd dense_hess(hess);

        auto f = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_f;
            Eigen::SparseMatrix<double> hess_f;
            compute_friction_term(
                x, constraints, grad_f, hess_f,
                /*compute_grad=*/true, /*compute_hess=*/false);
            return grad_f;
        };
        Eigen::MatrixXd hess_approx;
        fd::finite_jacobian(sigma, f, hess_approx);
        if (!fd::compare_hessian(hess, hess_approx)) {
            spdlog::error(
                "finite hessian check failed for friction "
                "(hess_L_inf_norm={:g} diff_L_inf_norm={:g})",
                dense_hess.lpNorm<Eigen::Infinity>(),
                (hess_approx - dense_hess).lpNorm<Eigen::Infinity>());
        }

        ///////////////////////////////////////////////////////////////////////
        // Auto. diff. check
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(sigma.size());
        Diff::D2MatrixXd V_diff = m_assembler.world_vertices(
            this->dofs_to_poses(Diff::d2vars(0, sigma)));

        // Contact constraints to friction constraints
        ipc::FrictionConstraints friction_constraint_set;
        ipc::construct_friction_constraint_set(
            V0, m_assembler.m_edges, m_assembler.m_faces, constraints,
            barrier_activation_distance(), barrier_stiffness(),
            coefficient_friction, friction_constraint_set);

        Diff::DDouble2 f_diff = compute_friction_potential(
            V0, V_diff, m_assembler.m_edges, m_assembler.m_faces,
            friction_constraint_set, static_friction_speed_bound * timestep());

        if (std::isfinite(f_diff.getHessian().sum())
            && !fd::compare_hessian(dense_hess, f_diff.getHessian())) {
            spdlog::error(
                "autodiff hessian check failed for friction "
                "(hess_L_inf_norm={:g} diff_L_inf_norm={:g})",
                dense_hess.lpNorm<Eigen::Infinity>(),
                (f_diff.getHessian() - dense_hess).lpNorm<Eigen::Infinity>());
        } else if (!std::isfinite(f_diff.getHessian().sum())) {
            spdlog::warn("autodiff hessian failed for friction");
        }
    }
#endif

} // namespace opt
} // namespace ccd
