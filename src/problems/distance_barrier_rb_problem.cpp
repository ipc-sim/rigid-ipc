#include "distance_barrier_rb_problem.hpp"

#include <finitediff.hpp>
// IPC Toolkit
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/ipc.hpp>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

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
        , body_energy_integration_method(DEFAULT_BODY_ENERGY_INTEGRATION_METHOD)
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

        body_energy_integration_method =
            params["rigid_body_problem"]["time_stepper"]
                .get<BodyEnergyIntegrationMethod>();
        RigidBodyProblem::settings(params["rigid_body_problem"]);

        if (friction_iterations == 0) {
            spdlog::info(
                "Disabling friction because friction iterations is zero");
            coefficient_friction = 0; // This disables all friction computation
        } else if (friction_iterations < 0) {
            friction_iterations = Constants::MAXIMUM_FRICTION_ITERATIONS;
        }

        min_distance = compute_min_distance(starting_point());
        if (min_distance < 0) {
            spdlog::info("init_min_distance=N/A");
        } else {
            spdlog::info("init_min_distance={:.8e}", min_distance);
        }
    }

    nlohmann::json DistanceBarrierRBProblem::settings() const
    {
        nlohmann::json json = RigidBodyProblem::settings();
        json["friction_iterations"] = friction_iterations;
        json["static_friction_speed_bound"] = static_friction_speed_bound;
        json["time_stepper"] = body_energy_integration_method;
        return json;
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
        tbb::parallel_for(size_t(0), num_bodies(), [&](size_t i) {
            m_assembler.m_rbs[i].pose_prev = m_assembler.m_rbs[i].pose;
        });
        // Update the stored poses and inital value for the solver
        update_dof();

        // Reset m_had_collision which will be filled in by has_collisions().
        m_had_collisions = false;
        m_num_contacts = 0;

        // Disable barriers if solve_collision == false
        this->m_use_barriers = solve_collisions;

        // Solve constraints updates the constraints and takes the step
        update_constraints();
        opt_result = solve_constraints();
        _has_intersections = take_step(opt_result.x);
        had_collisions = m_had_collisions;
    }

    void DistanceBarrierRBProblem::update_constraints()
    {
        PROFILE_POINT("DistanceBarrierRBProblem::update_constraints");
        PROFILE_START();

        RigidBodyProblem::update_constraints();

        ipc::Constraints collision_constraints;
        m_constraint.construct_constraint_set(
            m_assembler, poses_t0, collision_constraints);

        Eigen::SparseMatrix<double> hess;
        compute_barrier_term(
            x0, collision_constraints, grad_barrier_t0, hess,
            /*compute_grad=*/true, /*compute_hess=*/false);

        update_friction_constraints(collision_constraints, poses_t0);

        PROFILE_END();
    }

    void DistanceBarrierRBProblem::update_friction_constraints(
        const ipc::Constraints& collision_constraints,
        const physics::Poses<double>& poses)
    {
        PROFILE_POINT("DistanceBarrierRBProblem::update_friction_constraints");
        PROFILE_START();

        // The fricition constraints are constant through out the entire
        // lagging iteration.
        if (coefficient_friction > 0) {
            // Contact constraints to friction constraints
            friction_constraints.clear();
            Eigen::MatrixXd V0 = m_assembler.world_vertices(poses);
            ipc::construct_friction_constraint_set(
                V0, edges(), faces(), collision_constraints,
                barrier_activation_distance(), barrier_stiffness(),
                coefficient_friction, friction_constraints);
        }

        PROFILE_END();
    }

    opt::OptimizationResults DistanceBarrierRBProblem::solve_constraints()
    {
        opt::OptimizationResults opt_result;
        opt_result.x = starting_point();
        double momentum_balance, eps_d = 1e-2 * world_bbox_diagonal();
        int i = 0;
        int total_newton_iterations = 0;
        do {
            opt_result = solver().solve(opt_result.x);
            total_newton_iterations += opt_result.num_iterations;
            if (!opt_result.success) {
                break;
            }

            physics::Poses<double> poses = this->dofs_to_poses(opt_result.x);

            ipc::Constraints collision_constraints;
            m_constraint.construct_constraint_set(
                m_assembler, poses, collision_constraints);
            update_friction_constraints(collision_constraints, poses);

            Eigen::VectorXd grad_Ex, grad_Bx, grad_Dx;
            compute_energy_term(opt_result.x, grad_Ex);
            compute_barrier_term(opt_result.x, collision_constraints, grad_Bx);
            compute_friction_term(opt_result.x, grad_Dx);

            Eigen::VectorXd tmp =
                grad_Ex + barrier_stiffness() * grad_Bx + grad_Dx;
            tmp = is_dof_fixed().select(0, tmp);

            momentum_balance = tmp.norm();

            i++;

            spdlog::info(
                "friction_solve lagging_iteration={:d} momentum_balance={:g} "
                "eps_d={:g}",
                i, momentum_balance, eps_d);
        } while ((friction_iterations < 0 || i < friction_iterations)
                 && momentum_balance > eps_d);

        if (opt_result.success) {
            spdlog::info(
                "Finished friction solve after {:d} lagging iteration(s) and a "
                "momentum balance error of {:g}",
                i, momentum_balance);
        } else {
            spdlog::error(
                "Ending friction solve early because newton solve {:d} failed!",
                i);
        }

        opt_result.num_iterations = total_newton_iterations;
        return opt_result;
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
        double Ex =
            compute_energy_term(x, grad, hess, compute_grad, compute_hess);

        // The following is used to disable constraints if desired
        // (useful for testing).
        if (!m_use_barriers) {
            return Ex;
        }

        // Compute a common constraint set to use for contacts and friction
        // Start by updating the constraint set
        ipc::Constraints constraints;
        m_constraint.construct_constraint_set(
            m_assembler, this->dofs_to_poses(x), constraints);

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
        double Bx = compute_barrier_term(
            x, constraints, grad_Bx, hess_Bx, compute_grad, compute_hess);

        // D(x) is the friction potential (Equation 15 in the IPC paper)
        Eigen::VectorXd grad_Dx;
        Eigen::SparseMatrix<double> hess_Dx;
        double Dx = compute_friction_term(
            x, grad_Dx, hess_Dx, compute_grad, compute_hess);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                // Energy gradient is checked in compute_energy_term()
                check_grad_barrier(x, constraints, grad_Bx);
                // Friction gradient is checked in compute_friction_term()
            }
            if (compute_hess) {
                // Energy hessian is checked in compute_energy_term()
                check_hess_barrier(x, constraints, hess_Bx);
                // Friction hessian is checked in compute_friction_term()
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
        PROFILE_POINT("DistanceBarrierRBProblem::compute_energy_term");
        PROFILE_START();

        typedef AutodiffType<Eigen::Dynamic, /*maxN=*/6> Diff;

        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        int pos_ndof = physics::Pose<double>::dim_to_pos_ndof(dim());
        int rot_ndof = physics::Pose<double>::dim_to_rot_ndof(dim());

        Eigen::VectorXd energies(num_bodies());
        if (compute_grad) {
            grad.resize(x.size());
        }
        tbb::concurrent_vector<Eigen::Triplet<double>> hess_triplets;
        if (compute_hess) {
            // Hessian is a block diagonal with (ndof x ndof) blocks
            hess_triplets.reserve(num_bodies() * ndof * ndof);
        }

        const std::vector<physics::Pose<double>> poses = this->dofs_to_poses(x);
        assert(poses.size() == num_bodies());

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

                Diff::DDouble2 dExi = compute_body_energy<Diff::DDouble2>(
                    body, pose_diff, grad_barrier_t0.segment(i * ndof, ndof));

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

                Diff::DDouble1 dExi = compute_body_energy<Diff::DDouble1>(
                    body, pose_diff, grad_barrier_t0.segment(i * ndof, ndof));

                energies[i] = dExi.getValue();
                gradi = dExi.getGradient();
            } else {
                energies[i] = compute_body_energy<double>(
                    body, pose, grad_barrier_t0.segment(i * ndof, ndof));
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

        PROFILE_END();

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
        const physics::RigidBody& body,
        const physics::Pose<T>& pose,
        const Eigen::VectorX6d& grad_barrier_t0)
    {
        // NOTE: t0 suffix indicates the current value not the inital value
        double h = timestep();

        T energy(0.0);

        // Linear energy
        if (!body.is_dof_fixed.head(pose.pos_ndof()).all()) {
            Eigen::VectorX3<T> q = pose.position;
            const Eigen::VectorX3d& q_t0 = body.pose.position;
            const Eigen::VectorX3d& qdot_t0 = body.velocity.position;

            Eigen::VectorX3d a_t0;
            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                a_t0.setZero(pose.pos_ndof());
                break;
            case IMPLICIT_NEWMARK:
                a_t0 = 0.5 * grad_barrier_t0.head(pose.pos_ndof()) / body.mass;
                break;
            }
            Eigen::VectorX3d qdotdot_t0 =
                gravity + body.force.position / body.mass + a_t0;

            // ½mqᵀq - mqᵀ(qᵗ + h(q̇ᵗ + h(g + f/m + ½∇B(qᵗ)/m)))
            energy += 0.5 * body.mass * q.dot(q)
                - body.mass * q.dot((q_t0 + h * (qdot_t0 + h * (qdotdot_t0))));
        }

        // Rotational energy
        if (!body.is_dof_fixed.tail(pose.rot_ndof()).all()) {
            if (dim() == 3) {
                Eigen::Matrix3<T> Q = pose.construct_rotation_matrix();
                Eigen::Matrix3d Q_t0 = body.pose.construct_rotation_matrix();

                // Eigen::Matrix3d Qdot_t0 = Q_t0 *
                // Eigen::Hat(body.velocity.rotation);
                Eigen::Matrix3d Qdot_t0 = body.Qdot;

                const Eigen::VectorX3d& I = body.moment_of_inertia;
                Eigen::DiagonalMatrix<double, 3> J(
                    0.5 * (-I.x() + I.y() + I.z()), //
                    0.5 * (I.x() - I.y() + I.z()),  //
                    0.5 * (I.x() + I.y() - I.z()));

                Eigen::Matrix3d A_t0;
                switch (body_energy_integration_method) {
                case IMPLICIT_EULER:
                    A_t0.setZero();
                    break;
                case IMPLICIT_NEWMARK: {
                    Eigen::DiagonalMatrix<double, 3> J_inv(
                        2 / (-I.x() + I.y() + I.z()),
                        2 / (I.x() - I.y() + I.z()),
                        2 / (I.x() + I.y() - I.z()));
                    // A_t0 = 0.5 * ...;
                    throw NotImplementedError(
                        "Implicit Newmark not implemented in 3D!");
                    break;
                }
                }

                // ½tr(QJQᵀ) - tr(Q(J(Qᵗ + hQ̇ᵗ + h²Aᵗ)ᵀ + h²[τ]))
                Eigen::Matrix3d Tau =
                    Q_t0.transpose() * Eigen::Hat(body.force.rotation);
                energy += 0.5 * (Q * J * Q.transpose()).trace();
                energy -=
                    (Q * J * (Q_t0 + h * (Qdot_t0 + h * A_t0)).transpose())
                        .trace();
                energy += h * h * (Q * Tau).trace();
            } else {
                assert(pose.rot_ndof() == 1);
                T theta = pose.rotation[0];
                double theta_t0 = body.pose.rotation[0];
                double theta_dot_t0 = body.velocity.rotation[0];
                double tau = body.force.rotation[0];

                double I = body.moment_of_inertia[0];

                double a_t0;
                switch (body_energy_integration_method) {
                case IMPLICIT_EULER:
                    a_t0 = 0;
                    break;
                case IMPLICIT_NEWMARK:
                    a_t0 = 0.5 * grad_barrier_t0.tail(pose.rot_ndof())[0] / I;
                    break;
                }
                // θ̈ = τ/I + ½∇B(θᵗ)/I
                double theta_ddot_t0 = tau / I + a_t0;

                // ½Iθ² - Iθ(θᵗ + h(θ̇ᵗ + hθ̈ᵗ))
                double theta_hat =
                    theta_t0 + h * (theta_dot_t0 + h * theta_ddot_t0);
                energy += 0.5 * I * theta * theta - I * theta * theta_hat;
            }
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
        physics::Poses<double> poses = this->dofs_to_poses(x);
        ipc::Constraints constraints;
        m_constraint.construct_constraint_set(m_assembler, poses, constraints);
        num_constraints = constraints.num_constraints();

        m_num_contacts = std::max(m_num_contacts, num_constraints);

        spdlog::debug(
            "problem={} num_vertex_vertex_constraint={:d} "
            "num_edge_vertex_constraints={:d} num_edge_edge_constraints={:d} "
            "num_face_vertex_constraints={:d}",
            name(), constraints.vv_constraints.size(),
            constraints.ev_constraints.size(),
            constraints.ee_constraints.size(),
            constraints.fv_constraints.size());

        double Bx = compute_barrier_term(
            x, constraints, grad, hess, compute_grad, compute_hess);

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

    // NOTE: Avoid duplicate profile points by declaring them outside the
    // templated function.
    NAMED_PROFILE_POINT(
        "DistanceBarrierRBProblem::add_constraint_barrier:hessian",
        ADD_CONSTRAINT_BARRIER_HESS);
    NAMED_PROFILE_POINT(
        "DistanceBarrierRBProblem::add_constraint_barrier:gradient",
        ADD_CONSTRAINT_BARRIER_GRAD);
    NAMED_PROFILE_POINT(
        "DistanceBarrierRBProblem::add_constraint_barrier:value",
        ADD_CONSTRAINT_BARRIER_VAL);
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
            PROFILE_START(ADD_CONSTRAINT_BARRIER_HESS);

            Diff::DDouble2 dBxi = distance_barrier<Diff::DDouble2>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();
            Eigen::MatrixXd hessi = dBxi.getHessian();
            // Project dense block to make assembled matrix PSD
            hessi = Eigen::project_to_psd(hessi);
            // Add global triplets of the local values
            local_hessian_to_global_triplets(
                hessi, body_ids, ndof, hess_triplets);

            PROFILE_END(ADD_CONSTRAINT_BARRIER_HESS);
        } else if (compute_grad) {
            PROFILE_START(ADD_CONSTRAINT_BARRIER_GRAD);

            Diff::DDouble1 dBxi = distance_barrier<Diff::DDouble1>(sigma, rbc);
            Bx += dBxi.getValue();
            gradi = dBxi.getGradient();

            PROFILE_END(ADD_CONSTRAINT_BARRIER_GRAD);
        } else {
            PROFILE_START(ADD_CONSTRAINT_BARRIER_VAL);
            Bx += distance_barrier<double>(sigma, rbc);
            PROFILE_END(ADD_CONSTRAINT_BARRIER_VAL);
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
        PROFILE_POINT("DistanceBarrierRBProblem::compute_barrier_term");
        PROFILE_START();

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

        PROFILE_END();
        return Bx;
    }

    // Compute the derivatives of the function that maps from rigid body DOF
    // to vertex positions.
    Eigen::MatrixXd DistanceBarrierRBProblem::rigid_dof_to_vertices(
        const Eigen::VectorXd& x,
        Eigen::MatrixXd& jac,
        std::vector<Eigen::SparseMatrix<double>>& hess,
        bool compute_jac,
        bool compute_hess)
    {
        PROFILE_POINT("DistanceBarrierRBProblem::rigid_dof_to_vertices");
        PROFILE_START();

        // V: Rᵐ ↦ Rⁿ (vertices flattened rowwise)
        int m = x.size();
        int n = num_vertices() * dim();
        int rb_ndof = physics::Pose<double>::dim_to_ndof(dim());
        int rb_pos_ndof = physics::Pose<double>::dim_to_pos_ndof(dim());
        int rb_rot_ndof = physics::Pose<double>::dim_to_rot_ndof(dim());

        // We will only use auto diff to compute the derivatives of the rotation
        // matrix.
        typedef AutodiffType<Eigen::Dynamic, /*maxN=*/3> Diff;

        if (compute_jac) {
            // ∇ V(x): Rᵐ ↦ Rⁿˣᵐ
            jac = Eigen::MatrixXd::Zero(n, m);
        }
        if (compute_hess) {
            // ∇²V(x): Rᵐ ↦ Rⁿˣᵐˣᵐ
            hess.resize(n, Eigen::SparseMatrix<double>(m, m));
        }

        Eigen::MatrixXd V(num_vertices(), dim());
        tbb::parallel_for(size_t(0), num_bodies(), [&](size_t rb_i) {
            // Activate autodiff with the correct number of variables.
            Diff::activate(rb_rot_ndof);

            const physics::RigidBody& rb = m_assembler.m_rbs[rb_i];
            // Index of ribid bodies first vertex in the global vertices
            long rb_v0_i = m_assembler.m_body_vertex_id[rb_i];

            const auto& p = x.segment(rb_i * rb_ndof, rb_pos_ndof);
            const auto& r =
                x.segment(rb_i * rb_ndof + rb_pos_ndof, rb_rot_ndof);

            if (compute_hess) {
                auto R = construct_rotation_matrix(
                    Eigen::VectorX3<Diff::DDouble2>(Diff::d2vars(0, r)));
                Diff::D2MatrixXd V_diff = rb.vertices * R.transpose();

                for (int i = 0; i < V_diff.rows(); i++) {
                    for (int j = 0; j < V_diff.cols(); j++) {
                        V(rb_v0_i + i, j) = V_diff(i, j).getValue() + p(j);

                        // Fill in gradient of V(i, j) (∈ R⁶ for 3D)
                        int vij_flat = (rb_v0_i + i) * V_diff.cols() + j;
                        jac(vij_flat, rb_i * rb_ndof + j) = 1;
                        jac.block(
                            /*i=*/vij_flat, /*j=*/rb_i * rb_ndof + rb_pos_ndof,
                            /*p=*/1, /*q=*/rb_rot_ndof) =
                            V_diff(i, j).getGradient().transpose();

                        // Fill in hessian of V(i, j) (∈ R⁶ˣ⁶ for 3D)
                        // Hessian of position is zero
                        std::vector<Eigen::Triplet<double>> hess_triplets;
                        const auto& local_hess = V_diff(i, j).getHessian();
                        hess_triplets.reserve(rb_rot_ndof * rb_rot_ndof);
                        for (int k = 0; k < rb_rot_ndof; k++) {
                            for (int l = 0; l < rb_rot_ndof; l++) {
                                hess_triplets.emplace_back(
                                    rb_i * rb_ndof + rb_pos_ndof + k,
                                    rb_i * rb_ndof + rb_pos_ndof + l,
                                    local_hess(k, l));
                            }
                        }
                        hess[vij_flat].setFromTriplets(
                            hess_triplets.begin(), hess_triplets.end());
                    }
                }
            } else if (compute_jac) {
                auto R = construct_rotation_matrix(
                    Eigen::VectorX3<Diff::DDouble1>(Diff::d1vars(0, r)));
                Diff::D1MatrixXd V_diff = rb.vertices * R.transpose();

                for (int i = 0; i < V_diff.rows(); i++) {
                    for (int j = 0; j < V_diff.cols(); j++) {
                        V(rb_v0_i + i, j) = V_diff(i, j).getValue() + p(j);

                        // Fill in gradient of V(i, j) (∈ R⁶ for 3D)
                        int vij_flat = (rb_v0_i + i) * V_diff.cols() + j;
                        jac(vij_flat, rb_i * rb_ndof + j) = 1;
                        jac.block(
                            /*i=*/vij_flat, /*j=*/rb_i * rb_ndof + rb_pos_ndof,
                            /*p=*/1, /*q=*/rb_rot_ndof) =
                            V_diff(i, j).getGradient().transpose();
                    }
                }
            } else {
                V.block(rb_v0_i, 0, rb.vertices.rows(), rb.dim()) =
                    rb.world_vertices(physics::Pose<double>(
                        x.segment(rb_i * rb_ndof, rb_ndof)));
            }
        });

        assert(
            (V - m_assembler.world_vertices(this->dofs_to_poses(x))).norm()
            < 1e-12);

        PROFILE_END();
        return V;
    }

    double DistanceBarrierRBProblem::compute_friction_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        if (coefficient_friction <= 0) {
            grad.setZero(x.size());
            hess = Eigen::SparseMatrix<double>(x.size(), x.size());
            return 0;
        }

        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term",
            COMPUTE_FRICTION_TERM);
        NAMED_PROFILE_POINT("compute_friction_potential", COMPUTE_POTENTIAL);
        NAMED_PROFILE_POINT(
            "compute_friction_potential_gradient", COMPUTE_POTENTIAL_GRAD);
        NAMED_PROFILE_POINT(
            "compute_friction_potential_hessian", COMPUTE_POTENTIAL_HESS);

        PROFILE_START(COMPUTE_FRICTION_TERM);

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);

        // Compute V(x)
        Eigen::MatrixXd jac_V;
        std::vector<Eigen::SparseMatrix<double>> hess_V;
        Eigen::MatrixXd V = rigid_dof_to_vertices(
            x, jac_V, hess_V, compute_grad || compute_hess, compute_hess);

        // Compute the surface friction potential
        PROFILE_START(COMPUTE_POTENTIAL);
        double friction_potential = compute_friction_potential(
            V0, V, edges(), faces(), friction_constraints,
            static_friction_speed_bound * timestep());
        PROFILE_END(COMPUTE_POTENTIAL);

        Eigen::VectorXd grad_f;
        if (compute_grad || compute_hess) {
            PROFILE_START(COMPUTE_POTENTIAL_GRAD);
            // The gradient is also needed for the full hessian
            grad_f = compute_friction_potential_gradient(
                V0, V, edges(), faces(), friction_constraints,
                static_friction_speed_bound * timestep());
            PROFILE_END(COMPUTE_POTENTIAL_GRAD);
        }

        Eigen::SparseMatrix<double> hess_f;
        if (compute_hess) {
            PROFILE_START(COMPUTE_POTENTIAL_HESS);
            hess_f = compute_friction_potential_hessian(
                V0, V, edges(), faces(), friction_constraints,
                static_friction_speed_bound * timestep(),
                /*project_to_psd=*/false);
            PROFILE_END(COMPUTE_POTENTIAL_HESS);
        }

        // Apply the chain rule
        if (compute_grad) {
            // ∇ₓf(g(x)) = ∇ᵤf(u=g(x)) * ∇ₓV(x)
            grad = jac_V.transpose() * grad_f; // grad = [∇ₓf(V(x))]ᵀ
        }
        if (compute_hess) {
            Eigen::MatrixXd dense_hess = jac_V.transpose() * hess_f * jac_V;
            for (int i = 0; i < hess_V.size(); i++) {
                dense_hess += hess_V[i] * grad_f[i];
            }
            hess = Eigen::project_to_psd(dense_hess).sparseView();
        }

        PROFILE_END(COMPUTE_FRICTION_TERM);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_grad_friction(x, grad);
            }
            if (compute_hess) {
                check_hess_friction(x, hess);
            }
            is_checking_derivative = false;
        }
#endif

        return friction_potential;
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
            d_edge0_vertex0, d_edge0_vertex1, d_edge1_vertex0, d_edge1_vertex1);
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
                barrier_activation_distance(), d.getValue());
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
                barrier_activation_distance(), d.getValue());
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
        const Eigen::VectorXd& sigma, const Eigen::VectorXd& grad)
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
            return compute_friction_term(x);
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
        Diff::DDouble1 f_diff = compute_friction_potential(
            V0, V_diff, edges(), faces(), friction_constraints,
            static_friction_speed_bound * timestep());

        if (std::isfinite(f_diff.getGradient().sum())
            && !fd::compare_gradient(grad, f_diff.getGradient())) {
            spdlog::error("autodiff gradient check failed for friction");
        }
    }

    void DistanceBarrierRBProblem::check_hess_friction(
        const Eigen::VectorXd& sigma, const Eigen::SparseMatrix<double>& hess)
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
            compute_friction_term(x, grad_f);
            return grad_f;
        };
        Eigen::MatrixXd hess_approx;
        fd::finite_jacobian(sigma, f, hess_approx);
        hess_approx = Eigen::project_to_psd(hess_approx);
        if (!fd::compare_hessian(hess, hess_approx, 1e-2)) {
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

        Diff::DDouble2 f_diff = compute_friction_potential(
            V0, V_diff, edges(), faces(), friction_constraints,
            static_friction_speed_bound * timestep());

        Eigen::MatrixXd hess_autodiff =
            Eigen::project_to_psd(f_diff.getHessian());

        if (std::isfinite(hess_autodiff.sum())) {
            if (!fd::compare_hessian(dense_hess, hess_autodiff)) {
                spdlog::error(
                    "autodiff hessian check failed for friction "
                    "(hess_L_inf_norm={:g} diff_L_inf_norm={:g})",
                    dense_hess.lpNorm<Eigen::Infinity>(),
                    (hess_autodiff - dense_hess).lpNorm<Eigen::Infinity>());
            }
        } else {
            spdlog::warn("autodiff hessian failed for friction");
        }
    }
#endif

} // namespace opt
} // namespace ccd
