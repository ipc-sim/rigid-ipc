#include "distance_barrier_rb_problem.hpp"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/ipc.hpp>

#ifdef WITH_DERIVATIVE_CHECK
#include <finitediff.hpp>
#endif

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
            spdlog::info(
                "init_min_distance>d̂={:.8e}", barrier_activation_distance());
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

    Eigen::VectorXi DistanceBarrierRBProblem::free_dof() const
    {
        const Eigen::VectorXb& is_dof_fixed = this->is_dof_fixed();
        std::vector<int> free_dofs;
        free_dofs.reserve(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed[i] && !is_dof_satisfied[i]) {
                free_dofs.push_back(i);
            }
        }
        return Eigen::Map<Eigen::VectorXi>(free_dofs.data(), free_dofs.size());
    }

    ////////////////////////////////////////////////////////////
    // Rigid Body Problem

    void DistanceBarrierRBProblem::simulation_step(
        bool& had_collisions, bool& _has_intersections, bool solve_collisions)
    {
        // Advance the poses, but leave the current pose unchanged for now.
        for (size_t i = 0; i < num_bodies(); i++) {
            m_assembler[i].pose_prev = m_assembler[i].pose;
        }

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

        init_augmented_lagrangian();

        PROFILE_END();
    }

    void DistanceBarrierRBProblem::update_friction_constraints(
        const ipc::Constraints& collision_constraints,
        const physics::Poses<double>& poses)
    {
        if (coefficient_friction <= 0) {
            return;
        }

        PROFILE_POINT("DistanceBarrierRBProblem::update_friction_constraints");
        PROFILE_START();

        // The fricition constraints are constant through out the entire
        // lagging iteration.
        friction_constraints.clear();
        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses);
        ipc::construct_friction_constraint_set(
            V0, edges(), faces(), collision_constraints,
            barrier_activation_distance(), barrier_stiffness(),
            coefficient_friction, friction_constraints);

        PROFILE_END();
    }

    void DistanceBarrierRBProblem::init_augmented_lagrangian()
    {
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        int pos_ndof = physics::Pose<double>::dim_to_pos_ndof(dim());
        int rot_ndof = physics::Pose<double>::dim_to_rot_ndof(dim());

        augmented_lagrangian_penalty = 1e6;
        augmented_lagrangian_multiplier.setZero(
            ndof * m_assembler.num_kinematic_bodies());

        x_pred = x0;
        for (int i = 0; i < num_bodies(); i++) {
            x_pred.segment(ndof * i, pos_ndof) +=
                timestep() * m_assembler[i].velocity.position;
            x_pred.segment(ndof * i + pos_ndof, rot_ndof) +=
                timestep() * m_assembler[i].velocity.rotation;
        }

        is_dof_satisfied.setZero(x0.size());
        for (int i = 0; i < num_bodies(); i++) {
            if (m_assembler[i].type == physics::RigidBodyType::STATIC) {
                is_dof_satisfied.segment(ndof * i, ndof).setOnes();
            }
        }
    }

    double DistanceBarrierRBProblem::compute_convergence_criteria(
        const Eigen::VectorXd& x) const
    {
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        double diff_x_pred_x = 0, diff_x_pred_x0 = 0;
        for (const auto& k : m_assembler.m_kinematic_body_ids) {
            diff_x_pred_x +=
                (x_pred.segment(ndof * k, ndof) - x.segment(ndof * k, ndof))
                    .squaredNorm();
            diff_x_pred_x0 +=
                (x_pred.segment(ndof * k, ndof) - x0.segment(ndof * k, ndof))
                    .squaredNorm();
        }
        return 1 - sqrt(diff_x_pred_x / diff_x_pred_x0);
    }

    void DistanceBarrierRBProblem::update_augmented_lagrangian(
        const Eigen::VectorXd& x)
    {
        int ndof = physics::Pose<double>::dim_to_ndof(dim());

        double eta_AL = compute_convergence_criteria(x);
        if (eta_AL >= 0.999) {
            // Fix the kinematic DoF that have converged
            for (int i = 0; i < num_bodies(); i++) {
                if (m_assembler[i].type == physics::RigidBodyType::KINEMATIC) {
                    is_dof_satisfied.segment(ndof * i, ndof).setOnes();
                }
            }
            return;
        }

        if (eta_AL < 0.99 && augmented_lagrangian_penalty < 1e8) {
            augmented_lagrangian_penalty *= 2;
        } else {
            for (int ki = 0; ki < m_assembler.num_kinematic_bodies(); ki++) {
                const auto& k = m_assembler.m_kinematic_body_ids[ki];
                augmented_lagrangian_multiplier.segment(ki * ndof, ndof) -=
                    augmented_lagrangian_penalty * sqrt(m_assembler[k].mass)
                    * (x.segment(ndof * k, ndof)
                       - x_pred.segment(ndof * k, ndof));
            }
        }
        spdlog::info(
            "κ_AL={:g} ||λ_AL||∞={}", augmented_lagrangian_penalty,
            augmented_lagrangian_multiplier.lpNorm<Eigen::Infinity>());
    }

    bool DistanceBarrierRBProblem::are_equality_constraints_satisfied(
        const Eigen::VectorXd& x) const
    {
        if (m_assembler.num_kinematic_bodies()) {
            return compute_convergence_criteria(x) >= 0.999;
        }
        return true;
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

    bool DistanceBarrierRBProblem::take_step(const Eigen::VectorXd& x)
    {
        min_distance = compute_min_distance(x);
        if (min_distance < 0) {
            spdlog::info("final_step min_distance=N/A");
        } else {
            spdlog::info("final_step min_distance={:.8e}", min_distance);
        }

        return RigidBodyProblem::take_step(x);
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
        Ex /= average_mass();
        if (compute_grad) {
            grad /= average_mass();
        }
        if (compute_hess) {
            hess /= average_mass();
        }

        Eigen::VectorXd grad_AL;
        Eigen::SparseMatrix<double> hess_AL;
        double ALx = compute_augmented_lagrangian(
            x, grad_AL, hess_AL, compute_grad, compute_hess);
        Ex += ALx / average_mass();
        if (compute_grad) {
            grad += grad_AL / average_mass();
        }
        if (compute_hess) {
            hess += hess_AL / average_mass();
        }

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

        // Sum all the potentials
        double kappa_over_avg_mass = barrier_stiffness() / average_mass();
        if (compute_grad) {
            grad += kappa_over_avg_mass * grad_Bx + grad_Dx / average_mass();
        }
        if (compute_hess) {
            hess += kappa_over_avg_mass * hess_Bx + hess_Dx / average_mass();
        }

        return Ex + kappa_over_avg_mass * Bx + Dx / average_mass();
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

        Eigen::VectorXd energies = Eigen::VectorXd::Zero(num_bodies());
        if (compute_grad) {
            grad.setZero(x.size());
        }
        tbb::concurrent_vector<Eigen::Triplet<double>> hess_triplets;
        if (compute_hess) {
            // Hessian is a block diagonal with (ndof x ndof) blocks
            hess_triplets.reserve(num_bodies() * ndof * ndof);
        }

        const std::vector<physics::Pose<double>> poses = this->dofs_to_poses(x);
        assert(poses.size() == num_bodies());

        // tbb::parallel_for(size_t(0), poses.size(), [&](size_t i) {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), poses.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                // Activate autodiff with the correct number of variables
                Diff::activate(ndof);

                for (long i = range.begin(); i != range.end(); ++i) {
                    const physics::Pose<double>& pose = poses[i];
                    const physics::RigidBody& body = m_assembler[i];

                    // Do not compute the body energy for static and kinematic
                    // bodies
                    if (body.type != physics::RigidBodyType::DYNAMIC) {
                        continue;
                    }

                    Eigen::VectorX6d gradi;

                    if (compute_hess) {
                        // Initialize autodiff variables
                        physics::Pose<Diff::DDouble2> pose_diff(
                            Diff::d2vars(0, pose.dof()));

                        Diff::DDouble2 dExi =
                            compute_body_energy<Diff::DDouble2>(
                                body, pose_diff,
                                grad_barrier_t0.segment(i * ndof, ndof));

                        energies[i] = dExi.getValue();
                        gradi = dExi.getGradient();
                        Eigen::MatrixXX6d hessi = dExi.getHessian();

                        // Project dense block to make assembled matrix PSD
                        hessi.topLeftCorner(pos_ndof, pos_ndof) =
                            Eigen::project_to_psd(Eigen::MatrixXX3d(
                                hessi.topLeftCorner(pos_ndof, pos_ndof)));
                        // NOTE: This is handled with Tikhonov regularization
                        // instead hessi.bottomRightCorner(rot_ndof, rot_ndof) =
                        //     Eigen::project_to_pd(
                        //         hessi.bottomRightCorner(rot_ndof, rot_ndof));

                        // Add the local hessian as triplets in the global
                        // hessian
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
                            compute_body_energy<Diff::DDouble1>(
                                body, pose_diff,
                                grad_barrier_t0.segment(i * ndof, ndof));

                        energies[i] = dExi.getValue();
                        gradi = dExi.getGradient();
                    } else {
                        energies[i] = compute_body_energy<double>(
                            body, pose,
                            grad_barrier_t0.segment(i * ndof, ndof));
                    }

                    if (compute_grad) {
                        grad.segment(i * ndof, ndof) = gradi;
                    }
                }
            });

        if (compute_hess) {
            NAMED_PROFILE_POINT(
                "DistanceBarrierRBProblem::compute_energy_term:"
                "assemble_hessian",
                ASSEMBLE_ENERGY_HESS);
            PROFILE_START(ASSEMBLE_ENERGY_HESS);

            // ∇²E: Rⁿ ↦ Rⁿˣⁿ
            hess.resize(x.size(), x.size());
            hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

            PROFILE_END(ASSEMBLE_ENERGY_HESS);
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

    double DistanceBarrierRBProblem::compute_augmented_lagrangian(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        int ndof = physics::Pose<double>::dim_to_ndof(dim());
        int pos_ndof = physics::Pose<double>::dim_to_pos_ndof(dim());
        int rot_ndof = physics::Pose<double>::dim_to_rot_ndof(dim());

        double potential = 0;
        if (compute_grad) {
            grad.setZero(x.size());
        }
        std::vector<Eigen::Triplet<double>> hess_triplets;
        if (compute_hess) {
            hess.resize(x.size(), x.size());
            hess_triplets.reserve(m_assembler.num_kinematic_bodies() * ndof);
        }

        bool all_kinematic_dof_satisfied = true;
        for (const auto& k : m_assembler.m_kinematic_body_ids) {
            all_kinematic_dof_satisfied &=
                is_dof_satisfied.segment(ndof * k, ndof).all();
        }

        if (all_kinematic_dof_satisfied) {
            return potential;
        }

        PROFILE_POINT("DistanceBarrierRBProblem::compute_augmented_lagrangian");
        PROFILE_START();

        const auto& kappa = augmented_lagrangian_penalty;
        for (int ki = 0; ki < m_assembler.num_kinematic_bodies(); ki++) {
            const auto& k = m_assembler.m_kinematic_body_ids[ki];

            double m = m_assembler[k].mass;
            const auto& lambda =
                augmented_lagrangian_multiplier.segment(ki * ndof, ndof);

            const auto& x_k = x.segment(k * ndof, ndof);
            const auto& x_pred_k = x_pred.segment(k * ndof, ndof);

            potential += kappa / 2 * m * (x_k - x_pred_k).squaredNorm()
                - sqrt(m) * lambda.dot(x_k - x_pred_k);
            if (compute_grad) {
                grad.segment(k * ndof, ndof) =
                    kappa * m * (x_k - x_pred_k) - sqrt(m) * lambda;
            }
            if (compute_hess) {
                for (int i = 0; i < ndof; i++) {
                    hess_triplets.emplace_back(
                        ndof * k + i, ndof * k + i, kappa * m);
                }
            }
        }

        if (compute_hess) {
            NAMED_PROFILE_POINT(
                "DistanceBarrierRBProblem::compute_augmented_lagrangian:"
                "assemble_hessian",
                ASSEMBLE_AL_HESS);
            PROFILE_START(ASSEMBLE_AL_HESS);

            hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

            PROFILE_END(ASSEMBLE_AL_HESS);
        }

        PROFILE_END();

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_augmented_lagrangian_gradient(x, grad);
            }
            if (compute_hess) {
                check_augmented_lagrangian_hessian(x, hess);
            }
            is_checking_derivative = false;
        }
#endif

        return potential;
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

        return Bx;
    }

    // Convert from a local hessian to the triplets in the global hessian
    template <typename DerivedLocalGradient, typename Triplets>
    void local_gradient_to_global_triplets(
        const Eigen::MatrixBase<DerivedLocalGradient>& local_gradient,
        const std::array<long, 2>& body_ids,
        int ndof,
        Triplets& triplets)
    {
        assert(local_gradient.size() == 2 * ndof);
        for (int b_i = 0; b_i < body_ids.size(); b_i++) {
            for (int dof_i = 0; dof_i < ndof; dof_i++) {
                double v = local_gradient(ndof * b_i + dof_i);
                int r = ndof * body_ids[b_i] + dof_i;
                triplets.emplace_back(r, /*c=*/0, v);
            }
        }
    }

    template <typename DerivedLocalHessian, typename Triplets>
    void local_hessian_to_global_triplets(
        const Eigen::MatrixBase<DerivedLocalHessian>& local_hessian,
        const std::array<long, 2>& body_ids,
        int ndof,
        Triplets& triplets)
    {
        assert(local_hessian.rows() == 2 * ndof);
        assert(local_hessian.cols() == 2 * ndof);
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

    // Apply the chain rule of f(V(x)) given ∇ᵥf(V) and ∇ₓV(x)
    void apply_chain_rule(
        const Eigen::VectorX12d& grad_f,
        const Eigen::MatrixXd& jac_V,
        const Eigen::MatrixXX12d& hess_f,
        const Eigen::MatrixXd& hess_V,
        const std::vector<long>& vertex_ids,
        const std::vector<uint8_t>& local_body_ids,
        const std::array<long, 2>& body_ids,
        const int dim,
        std::vector<Eigen::Triplet<double>>& grad_triplets,
        std::vector<Eigen::Triplet<double>>& hess_triplets,
        bool compute_grad,
        bool compute_hess)
    {
        if (!compute_grad && !compute_hess) {
            return;
        }

        const int rb_ndof = physics::Pose<double>::dim_to_ndof(dim);

        if (compute_grad) {
            // jac_Vi ∈ R^{4n × 2m}
            Eigen::VectorX12d grad = Eigen::VectorX12d::Zero(2 * rb_ndof);
            for (int i = 0; i < vertex_ids.size(); i++) {
                grad.segment(rb_ndof * local_body_ids[i], rb_ndof) +=
                    jac_V.middleRows(vertex_ids[i] * dim, dim).transpose()
                    * grad_f.segment(i * dim, dim);
            }

            local_gradient_to_global_triplets(
                grad, body_ids, rb_ndof, grad_triplets);
        }

        if (compute_hess) {
            // jac_Vi ∈ R^{4n × 2m}
            Eigen::MatrixXX12d jac_Vi =
                Eigen::MatrixXX12d::Zero(vertex_ids.size() * dim, 2 * rb_ndof);
            for (int i = 0; i < vertex_ids.size(); i++) {
                jac_Vi.block(
                    i * dim, local_body_ids[i] * rb_ndof, dim, rb_ndof) =
                    jac_V.middleRows(vertex_ids[i] * dim, dim);
            }

            // hess ∈ R^{2m × 2m}
            Eigen::MatrixXX12d hess = jac_Vi.transpose() * hess_f * jac_Vi;
            for (int i = 0; i < vertex_ids.size(); i++) {
                for (int j = 0; j < dim; j++) {
                    // Off diagaonal blocks are all zero because the derivative
                    // of a vertex of body A with body B is zero.
                    hess.block(
                        local_body_ids[i] * rb_ndof,
                        local_body_ids[i] * rb_ndof, rb_ndof, rb_ndof) +=
                        hess_V.middleRows(
                            rb_ndof * (vertex_ids[i] * dim + j), rb_ndof)
                        * grad_f[i * dim + j];
                }
            }

            hess = project_to_psd(hess);

            local_hessian_to_global_triplets(
                hess, body_ids, rb_ndof, hess_triplets);
        }
    }

    double DistanceBarrierRBProblem::compute_barrier_term(
        const Eigen::VectorXd& x,
        const ipc::Constraints& constraints,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        if (constraints.size() == 0) {
            grad.setZero(x.size());
            hess = Eigen::SparseMatrix<double>(x.size(), x.size());
            return 0;
        }

        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_barrier_term",
            COMPUTE_BARRIER_TERM);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_barrier_term:gradient",
            COMPUTE_BARRIER_GRAD);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_barrier_term:assemble_gradient",
            ASSEMBLE_BARRIER_GRAD);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_barrier_term:hessian",
            COMPUTE_BARRIER_HESS);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_barrier_term:assemble_hessian",
            ASSEMBLE_BARRIER_HESS);

        PROFILE_START(COMPUTE_BARRIER_TERM);

        int rb_ndof = physics::Pose<double>::dim_to_ndof(dim());
        if (compute_grad) {
            PROFILE_START(COMPUTE_BARRIER_GRAD);
            grad.setZero(x.size());
        }
        if (compute_hess) {
            PROFILE_START(COMPUTE_BARRIER_HESS);
            hess.resize(x.size(), x.size());
        }

        // Compute V(x)
        Eigen::MatrixXd jac_V, hess_V;
        Eigen::MatrixXd V = m_assembler.world_vertices_diff(
            x, jac_V, hess_V, compute_grad || compute_hess, compute_hess);

        double dhat = barrier_activation_distance();

        Eigen::VectorXd potentials(constraints.size());
        struct ThreadStorage {
            // double potential = 0;
            std::vector<Eigen::Triplet<double>> grad_triplets;
            std::vector<Eigen::Triplet<double>> hess_triplets;
        };
        typedef tbb::enumerable_thread_specific<ThreadStorage> LocalStorage;
        LocalStorage storages;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), constraints.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                // Get references to the local derivative storage
                LocalStorage::reference local_storage = storages.local();
                // double& potential = local_storage.potential;
                auto& grad_triplets = local_storage.grad_triplets;
                auto& hess_triplets = local_storage.hess_triplets;

                // potential = 0;
                if (compute_grad) {
                    grad_triplets.reserve(2 * rb_ndof * range.size());
                }
                if (compute_hess) {
                    hess_triplets.reserve(4 * rb_ndof * rb_ndof * range.size());
                }

                for (size_t ci = range.begin(); ci != range.end(); ++ci) {
                    const auto& constraint = constraints[ci];

                    // potential +=
                    potentials[ci] =
                        constraint.compute_potential(V, edges(), faces(), dhat);

                    Eigen::VectorX12d grad_B;
                    if (compute_grad || compute_hess) {
                        grad_B = constraint.compute_potential_gradient(
                            V, edges(), faces(), dhat);
                    }

                    Eigen::MatrixXX12d hess_B;
                    if (compute_hess) {
                        hess_B = constraint.compute_potential_hessian(
                            V, edges(), faces(), dhat,
                            /*project_to_psd=*/false);
                    }

                    apply_chain_rule(
                        grad_B, jac_V, hess_B, hess_V,
                        constraint.vertex_indices(edges(), faces()),
                        vertex_local_body_ids(constraints, ci),
                        body_ids(m_assembler, constraints, ci), dim(),
                        grad_triplets, hess_triplets, compute_grad,
                        compute_hess);
                }
            });

        // double potential = 0;
        double potential = potentials.sum();
        for (const auto& local_storage : storages) {
            // potential += local_storage.potential;
            if (compute_grad) {
                PROFILE_START(ASSEMBLE_BARRIER_GRAD);

                Eigen::SparseMatrix<double> tmp_grad(x.size(), 1);
                tmp_grad.setFromTriplets(
                    local_storage.grad_triplets.begin(),
                    local_storage.grad_triplets.end());
                grad += tmp_grad;

                PROFILE_END(ASSEMBLE_BARRIER_GRAD);
            }
            if (compute_hess) {
                PROFILE_START(ASSEMBLE_BARRIER_HESS);

                Eigen::SparseMatrix<double> tmp_hess(x.size(), x.size());
                tmp_hess.setFromTriplets(
                    local_storage.hess_triplets.begin(),
                    local_storage.hess_triplets.end());
                hess += tmp_hess;

                PROFILE_END(ASSEMBLE_BARRIER_HESS);
            }
        }

        // if not active, nothing happens
        PROFILE_END(COMPUTE_BARRIER_GRAD);
        PROFILE_END(COMPUTE_BARRIER_HESS);
        PROFILE_END(COMPUTE_BARRIER_TERM);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_barrier_gradient(x, constraints, grad);
            }
            if (compute_hess) {
                check_barrier_hessian(x, constraints, hess);
            }
            is_checking_derivative = false;
        }
#endif

        return potential;
    }

    template <typename RigidBodyConstraint, typename FrictionConstraint>
    double DistanceBarrierRBProblem::compute_friction_potential(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXd& jac_V,
        const Eigen::MatrixXd& hess_V,
        const FrictionConstraint& constraint,
        tbb::concurrent_vector<Eigen::Triplet<double>>& grad_triplets,
        tbb::concurrent_vector<Eigen::Triplet<double>>& hess_triplets,
        bool compute_grad,
        bool compute_hess)
    {
        // for each constraint:
        //     let m ∈ {6}, n ∈ {2, 3, 4}
        //     compute    V(x) ∈ R^{3n},
        //             ∇ₓ V(x) ∈ R^{3n × 2m},
        //             ∇ₓ²V(x) ∈ R^{3n × 2m × 2m}
        //     compute    D(V) ∈ R,
        //             ∇ᵥ D(V) ∈ R^{3n},
        //             ∇ᵥ²D(V) ∈ R^{3n × 3n}
        //     ∇ₓD(V(x)) = ∇ₓV(x)ᵀ∇ᵥD(V) ∈ R^{12x12}
        //     ∇ₓ²D(V(x)) = ∇ₓV(x)ᵀ∇ᵥ²D(V)∇ₓV(x) + ∑ᵢ ∇ₓᵢ²V * ∇ᵥD[i]
        //     local_to_global(∇ₓD(V(x)))
        //     local_to_global(project_to_psd(∇ₓ²D(V(x))))

        int rb_ndof = physics::Pose<double>::dim_to_ndof(dim());

        Eigen::VectorX12d grad_D;
        Eigen::MatrixXX12d hess_D;
        double epsv_times_h = static_friction_speed_bound * timestep();
        double Dx =
            constraint.compute_potential(U, edges(), faces(), epsv_times_h);
        if (compute_grad || compute_hess) {
            grad_D = constraint.compute_potential_gradient(
                U, edges(), faces(), epsv_times_h);
        }
        if (compute_hess) {
            hess_D = constraint.compute_potential_hessian(
                U, edges(), faces(), epsv_times_h, /*project_to_psd=*/false);
        }
        std::vector<long> vertex_ids =
            constraint.vertex_indices(edges(), faces());

        RigidBodyConstraint rbc(m_assembler, constraint);
        auto local_body_ids = rbc.vertex_local_body_ids();

        if (compute_grad) {
            // jac_Vi ∈ R^{4n × 2m}
            Eigen::VectorX12d grad = Eigen::VectorX12d::Zero(2 * rb_ndof);
            for (int i = 0; i < vertex_ids.size(); i++) {
                grad.segment(rb_ndof * local_body_ids[i], rb_ndof) +=
                    jac_V.middleRows(vertex_ids[i] * dim(), dim()).transpose()
                    * grad_D.segment(i * dim(), dim());
            }

            local_gradient_to_global_triplets(
                grad, rbc.body_ids(), rb_ndof, grad_triplets);
        }

        if (compute_hess) {
            // jac_Vi ∈ R^{4n × 2m}
            Eigen::MatrixXX12d jac_Vi =
                Eigen::MatrixXd::Zero(vertex_ids.size() * dim(), 2 * rb_ndof);
            for (int i = 0; i < vertex_ids.size(); i++) {
                jac_Vi.block(
                    i * dim(), local_body_ids[i] * rb_ndof, dim(), rb_ndof) =
                    jac_V.middleRows(vertex_ids[i] * dim(), dim());
            }

            // hess ∈ R^{2m × 2m}
            Eigen::MatrixXX12d hess = jac_Vi.transpose() * hess_D * jac_Vi;
            for (int i = 0; i < vertex_ids.size(); i++) {
                for (int j = 0; j < dim(); j++) {
                    // Off diagaonal blocks are all zero because the derivative
                    // of a vertex of body A with body B is zero.
                    hess.block(
                        local_body_ids[i] * rb_ndof,
                        local_body_ids[i] * rb_ndof, rb_ndof, rb_ndof) +=
                        hess_V.middleRows(
                            rb_ndof * (vertex_ids[i] * dim() + j), rb_ndof)
                        * grad_D[i * dim() + j];
                }
            }

            hess = project_to_psd(hess);

            local_hessian_to_global_triplets(
                hess, rbc.body_ids(), rb_ndof, hess_triplets);
        }

        return Dx;
    }

    double DistanceBarrierRBProblem::compute_friction_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess)
    {
        if (coefficient_friction <= 0 || friction_constraints.size() == 0) {
            grad.setZero(x.size());
            hess = Eigen::SparseMatrix<double>(x.size(), x.size());
            return 0;
        }

        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term",
            COMPUTE_FRICTION_TERM);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term:gradient",
            COMPUTE_FRICTION_GRAD);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term:assemble_gradient",
            ASSEMBLE_FRICTION_GRAD);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term:hessian",
            COMPUTE_FRICTION_HESS);
        NAMED_PROFILE_POINT(
            "DistanceBarrierRBProblem::compute_friction_term:assemble_hessian",
            ASSEMBLE_FRICTION_HESS);

        PROFILE_START(COMPUTE_FRICTION_TERM);

        tbb::concurrent_vector<Eigen::Triplet<double>> grad_triplets;
        tbb::concurrent_vector<Eigen::Triplet<double>> hess_triplets;
        int rb_ndof = physics::Pose<double>::dim_to_ndof(dim());
        if (compute_grad) {
            PROFILE_START(COMPUTE_FRICTION_GRAD);
            grad_triplets.reserve(2 * rb_ndof * friction_constraints.size());
        }
        if (compute_hess) {
            PROFILE_START(COMPUTE_FRICTION_HESS);
            hess_triplets.reserve(
                4 * rb_ndof * rb_ndof * friction_constraints.size());
        }

        // Compute V(x)
        Eigen::MatrixXd jac_V, hess_V;
        Eigen::MatrixXd V1 = m_assembler.world_vertices_diff(
            x, jac_V, hess_V, compute_grad || compute_hess, compute_hess);

        // absolute linear dislacement of each point
        Eigen::MatrixXd U = V1 - m_assembler.world_vertices(poses_t0);

        Eigen::VectorXd friction_potential(friction_constraints.size());

        size_t start_ci = 0;
        tbb::parallel_for(
            size_t(0), friction_constraints.vv_constraints.size(),
            [&](size_t ci) {
                friction_potential[start_ci + ci] =
                    compute_friction_potential<RigidBodyVertexVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.vv_constraints[ci], grad_triplets,
                        hess_triplets, compute_grad, compute_hess);
            });

        start_ci += friction_constraints.vv_constraints.size();
        tbb::parallel_for(
            size_t(0), friction_constraints.ev_constraints.size(),
            [&](size_t ci) {
                friction_potential[start_ci + ci] =
                    compute_friction_potential<RigidBodyEdgeVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.ev_constraints[ci], grad_triplets,
                        hess_triplets, compute_grad, compute_hess);
            });

        start_ci += friction_constraints.ev_constraints.size();
        tbb::parallel_for(
            size_t(0), friction_constraints.ee_constraints.size(),
            [&](size_t ci) {
                friction_potential[start_ci + ci] =
                    compute_friction_potential<RigidBodyEdgeEdgeConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.ee_constraints[ci], grad_triplets,
                        hess_triplets, compute_grad, compute_hess);
            });

        start_ci += friction_constraints.ee_constraints.size();
        tbb::parallel_for(
            size_t(0), friction_constraints.fv_constraints.size(),
            [&](size_t ci) {
                friction_potential[start_ci + ci] =
                    compute_friction_potential<RigidBodyFaceVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.fv_constraints[ci], grad_triplets,
                        hess_triplets, compute_grad, compute_hess);
            });

        if (compute_grad) {
            PROFILE_START(ASSEMBLE_FRICTION_GRAD);

            Eigen::SparseMatrix<double> sparse_grad(x.size(), 1);
            sparse_grad.setFromTriplets(
                grad_triplets.begin(), grad_triplets.end());
            grad = Eigen::VectorXd(sparse_grad);

            PROFILE_END(ASSEMBLE_FRICTION_GRAD);
            PROFILE_END(COMPUTE_FRICTION_GRAD);
        }
        if (compute_hess) {
            PROFILE_START(ASSEMBLE_FRICTION_HESS);

            hess.resize(x.size(), x.size());
            hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

            PROFILE_END(ASSEMBLE_FRICTION_HESS);
            PROFILE_END(COMPUTE_FRICTION_HESS);
        }

        PROFILE_END(COMPUTE_FRICTION_TERM);

#ifdef WITH_DERIVATIVE_CHECK
        if (!is_checking_derivative) {
            is_checking_derivative = true;
            if (compute_grad) {
                check_friction_gradient(x, grad);
            }
            if (compute_hess) {
                check_friction_hessian(x, hess);
            }
            is_checking_derivative = false;
        }
#endif

        return friction_potential.sum();
    }

    ///////////////////////////////////////////////////////////////////////////

    double DistanceBarrierRBProblem::compute_min_distance(
        const Eigen::VectorXd& x) const
    {
        physics::Poses<double> poses = this->dofs_to_poses(x);
        double min_distance =
            m_constraint.compute_minimum_distance(m_assembler, poses);
        return std::isfinite(min_distance) ? min_distance : -1;
    }

    bool DistanceBarrierRBProblem::has_collisions(
        const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j)
    {
        physics::Poses<double> poses_i = this->dofs_to_poses(x_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(x_j);
        bool collisions =
            m_constraint.has_active_collisions(m_assembler, poses_i, poses_j);
        m_had_collisions |= collisions;
        // m_use_barriers := solve_collisions
        return m_use_barriers ? collisions : false;
    }

    double DistanceBarrierRBProblem::compute_earliest_toi(
        const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j)
    {
        // m_use_barriers := solve_collisions
        // If we are not solve collisions then just compute if there was a
        // collision.
        if (!m_use_barriers) {
            this->has_collisions(x_i, x_j); // will set m_had_collisions
            return std::numeric_limits<double>::infinity();
        }

        physics::Poses<double> poses_i = this->dofs_to_poses(x_i);
        physics::Poses<double> poses_j = this->dofs_to_poses(x_j);
        double earliest_toi =
            m_constraint.compute_earliest_toi(m_assembler, poses_i, poses_j);
        m_had_collisions |= earliest_toi <= 1;
        return earliest_toi;
    }

#ifdef WITH_DERIVATIVE_CHECK
    // The following functions are used exclusivly to check that the
    // gradient and hessian match a finite difference version.

    void DistanceBarrierRBProblem::check_barrier_gradient(
        const Eigen::VectorXd& x,
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
        fd::finite_gradient(x, b, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx, 1e-3)) {
            spdlog::error("finite gradient check failed for barrier");
        }
    }

    void DistanceBarrierRBProblem::check_barrier_hessian(
        const Eigen::VectorXd& x,
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

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        // WARNING: The following check does not work well because the
        // different projections to PSD can affect results.
        Eigen::MatrixXd dense_hess(hess);
        auto b = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_b;
            Eigen::SparseMatrix<double> hess_b;
            compute_barrier_term(
                x, constraints, grad_b, hess_b,
                /*compute_grad=*/true, /*compute_hess=*/false);
            return grad_b;
        };
        Eigen::MatrixXd hess_approx;
        fd::finite_jacobian(x, b, hess_approx);
        // hess_approx = Eigen::project_to_psd(hess_approx);
        if (!fd::compare_jacobian(
                hess, hess_approx, Constants::FINITE_DIFF_TEST)) {
            spdlog::error("finite hessian check failed for barrier");
        }
    }

    void DistanceBarrierRBProblem::check_friction_gradient(
        const Eigen::VectorXd& x, const Eigen::VectorXd& grad)
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
        fd::finite_gradient(x, f, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx)) {
            spdlog::error("finite gradient check failed for friction");
        }

        ///////////////////////////////////////////////////////////////////////
        // Auto. diff. check
        typedef AutodiffType<Eigen::Dynamic> Diff;
        Diff::activate(x.size());
        Diff::D1MatrixXd V_diff =
            m_assembler.world_vertices(this->dofs_to_poses(Diff::d1vars(0, x)));

        Eigen::MatrixXd V0 = m_assembler.world_vertices(poses_t0);
        Diff::DDouble1 f_diff = ipc::compute_friction_potential(
            V0, V_diff, edges(), faces(), friction_constraints,
            static_friction_speed_bound * timestep());

        if (std::isfinite(f_diff.getGradient().sum())
            && !fd::compare_gradient(grad, f_diff.getGradient())) {
            spdlog::error("autodiff gradient check failed for friction");
        }
    }

    void DistanceBarrierRBProblem::check_friction_hessian(
        const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& hess)
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
        Eigen::MatrixXd V1 = m_assembler.world_vertices(this->dofs_to_poses(x));
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
        fd::finite_jacobian(x, f, hess_approx);
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
        Diff::activate(x.size());
        Diff::D2MatrixXd V_diff =
            m_assembler.world_vertices(this->dofs_to_poses(Diff::d2vars(0, x)));

        Diff::DDouble2 f_diff = ipc::compute_friction_potential(
            V0, V_diff, edges(), faces(), friction_constraints,
            static_friction_speed_bound * timestep());

        Eigen::MatrixXd hess_autodiff =
            Eigen::project_to_psd(f_diff.getHessian());

        if (std::isfinite(hess_autodiff.sum())) {
            if (!fd::compare_hessian(dense_hess, hess_autodiff, 1e-3)) {
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

    void DistanceBarrierRBProblem::check_augmented_lagrangian_gradient(
        const Eigen::VectorXd& x, const Eigen::VectorXd& grad)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        for (int i = 0; i < grad.size(); i++) {
            if (!std::isfinite(grad(i))) {
                spdlog::error("augmented lagrangian gradient is not finite");
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        auto AL = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_AL;
            Eigen::SparseMatrix<double> hess_AL;
            return compute_augmented_lagrangian(
                x, grad_AL, hess_AL,
                /*compute_grad=*/false, /*compute_hess=*/false);
        };
        Eigen::VectorXd grad_approx;
        fd::finite_gradient(x, AL, grad_approx);
        if (!fd::compare_gradient(grad, grad_approx)) {
            spdlog::error(
                "finite gradient check failed for augmented lagrangian");
        }
    }

    void DistanceBarrierRBProblem::check_augmented_lagrangian_hessian(
        const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& hess)
    {
        ///////////////////////////////////////////////////////////////////////
        // Check that everything went well
        typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;
        for (int k = 0; k < hess.outerSize(); ++k) {
            for (Iterator it(hess, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    spdlog::error("augmented lagrangian hessian is not finite");
                    return;
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Finite difference check
        Eigen::MatrixXd dense_hess(hess);
        auto AL = [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd grad_AL;
            Eigen::SparseMatrix<double> hess_AL;
            compute_augmented_lagrangian(
                x, grad_AL, hess_AL,
                /*compute_grad=*/true, /*compute_hess=*/false);
            return grad_AL;
        };
        Eigen::MatrixXd hess_approx;
        fd::finite_jacobian(x, AL, hess_approx);
        if (!fd::compare_jacobian(hess, hess_approx)) {
            spdlog::error(
                "finite hessian check failed for augmented lagrangian");
        }
    }
#endif

} // namespace opt
} // namespace ccd
