#include "distance_barrier_rb_problem.hpp"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/ipc.hpp>

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
#include <finitediff.hpp>
#endif

#include <constants.hpp>
#include <geometry/distance.hpp>
#include <solvers/solver_factory.hpp>
#include <utils/not_implemented_error.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ipc::rigid {

DistanceBarrierRBProblem::DistanceBarrierRBProblem()
    : m_barrier_stiffness(1)
    , min_distance(-1)
    , m_had_collisions(false)
    , static_friction_speed_bound(1e-3)
    , friction_iterations(1)
    , body_energy_integration_method(DEFAULT_BODY_ENERGY_INTEGRATION_METHOD)
{
}

bool DistanceBarrierRBProblem::settings(const nlohmann::json& params)
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
    bool success = RigidBodyProblem::settings(params["rigid_body_problem"]);
    if (!success) {
        return false;
    }

    if (friction_iterations == 0) {
        spdlog::info("Disabling friction because friction iterations is zero");
        coefficient_friction = 0; // This disables all friction computation
    }
    // else if (friction_iterations < 0) {
    //     friction_iterations = Constants::MAXIMUM_FRICTION_ITERATIONS;
    // }

    min_distance = compute_min_distance(starting_point());
    if (min_distance < 0) {
        spdlog::info(
            "init_min_distance>d̂+dmin={:.8e}",
            barrier_activation_distance()
                + m_constraint.minimum_separation_distance);
    } else {
        spdlog::info("init_min_distance={:.8e}", min_distance);
    }

    return true;
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
    const VectorXb& is_dof_fixed = this->is_dof_fixed();
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
        m_assembler[i].velocity_prev = m_assembler[i].velocity;
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
    step_kinematic_bodies();
    had_collisions = m_had_collisions;
}

void DistanceBarrierRBProblem::update_constraints()
{
    PROFILE_POINT("DistanceBarrierRBProblem::update_constraints");
    PROFILE_START();

    RigidBodyProblem::update_constraints();

    Constraints collision_constraints;
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
    const Constraints& collision_constraints, const PosesD& poses)
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
    construct_friction_constraint_set(
        V0, edges(), faces(), collision_constraints,
        barrier_activation_distance(), barrier_stiffness(),
        coefficient_friction, friction_constraints);

    PROFILE_END();
}

void DistanceBarrierRBProblem::init_augmented_lagrangian()
{
    int ndof = PoseD::dim_to_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());
    int rot_ndof = PoseD::dim_to_rot_ndof(dim());

    linear_augmented_lagrangian_penalty = 1e3;
    angular_augmented_lagrangian_penalty = 1e3;
    size_t num_kinematic_bodies = m_assembler.count_kinematic_bodies();
    linear_augmented_lagrangian_multiplier.setZero(
        pos_ndof * num_kinematic_bodies);
    angular_augmented_lagrangian_multiplier.setZero(
        rot_ndof * num_kinematic_bodies, rot_ndof);

    for (int i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::KINEMATIC
            && m_assembler[i].kinematic_max_time < 0) {
            m_assembler[i].convert_to_static();
        }
    }

    x_pred = x0;
    for (int i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].kinematic_poses.size()) {
            const PoseD& pose = m_assembler[i].kinematic_poses.front();
            // Kinematic position
            x_pred.segment(ndof * i, pos_ndof) = pose.position;
            // Kinematic rotation
            x_pred.segment(ndof * i + pos_ndof, rot_ndof) = pose.rotation;
        } else {
            // Kinematic position
            x_pred.segment(ndof * i, pos_ndof) +=
                timestep() * m_assembler[i].velocity.position;
            // Kinematic rotation
            x_pred.segment(ndof * i + pos_ndof, rot_ndof) +=
                timestep() * m_assembler[i].velocity.rotation;
        }
    }

    is_dof_satisfied.setZero(x0.size());
    for (int i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::STATIC) {
            is_dof_satisfied.segment(ndof * i, ndof).setOnes();
        }
    }
}

void DistanceBarrierRBProblem::step_kinematic_bodies()
{
    for (int i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
            if (m_assembler[i].kinematic_max_time < 0) {
                m_assembler[i].convert_to_static();
            } else {
                m_assembler[i].kinematic_max_time -= timestep();
                if (m_assembler[i].kinematic_poses.size()) {
                    m_assembler[i].kinematic_poses.pop_front();
                }
            }
        }
    }
}

inline DiagonalMatrix3d compute_J(const VectorMax3d& I)
{
    return DiagonalMatrix3d(
        0.5 * (-I.x() + I.y() + I.z()), //
        0.5 * (I.x() - I.y() + I.z()),  //
        0.5 * (I.x() + I.y() - I.z()));
}

inline DiagonalMatrix3d compute_Jinv(const VectorMax3d& I)
{
    return DiagonalMatrix3d(
        2 / (-I.x() + I.y() + I.z()), //
        2 / (I.x() - I.y() + I.z()),  //
        2 / (I.x() + I.y() - I.z()));
}

inline DiagonalMatrix3d compute_Jsqrt(const VectorMax3d& I)
{
    assert(0.5 * (-I.x() + I.y() + I.z()) >= 0);
    assert(0.5 * (I.x() - I.y() + I.z()) >= 0);
    assert(0.5 * (I.x() + I.y() - I.z()) >= 0);
    return DiagonalMatrix3d(
        sqrt(std::max(0.5 * (-I.x() + I.y() + I.z()), 0.0)),
        sqrt(std::max(0.5 * (I.x() - I.y() + I.z()), 0.0)),
        sqrt(std::max(0.5 * (I.x() + I.y() - I.z()), 0.0)));
}

double DistanceBarrierRBProblem::compute_linear_augment_lagrangian_progress(
    const Eigen::VectorXd& x) const
{
    int ndof = PoseD::dim_to_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());

    double a = 0, b = 0;

    for (size_t i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
            a += (x_pred.segment(ndof * i, pos_ndof)
                  - x.segment(ndof * i, pos_ndof))
                     .squaredNorm();
            b += (x_pred.segment(ndof * i, pos_ndof)
                  - x0.segment(ndof * i, pos_ndof))
                     .squaredNorm();
        }
    }

    if (a == 0 && b == 0) {
        return 1;
    }
    return 1 - sqrt(a / (b == 0 ? 1 : b));
}

double DistanceBarrierRBProblem::compute_angular_augment_lagrangian_progress(
    const Eigen::VectorXd& x) const
{
    int ndof = PoseD::dim_to_ndof(dim());
    int rot_ndof = PoseD::dim_to_rot_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());

    double a = 0, b = 0;
    for (size_t i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
            size_t ri = ndof * i + pos_ndof;
            if (dim() == 2) {
                a += (x_pred.segment(ri, rot_ndof) - x.segment(ri, rot_ndof))
                         .squaredNorm();
                b += (x_pred.segment(ri, rot_ndof) - x0.segment(ri, rot_ndof))
                         .squaredNorm();
            } else {
                auto Q_pred = construct_rotation_matrix(
                    VectorMax3d(x_pred.segment(ri, rot_ndof)));
                auto Q = construct_rotation_matrix(
                    VectorMax3d(x.segment(ri, rot_ndof)));
                auto Q0 = construct_rotation_matrix(
                    VectorMax3d(x0.segment(ri, rot_ndof)));
                a += (Q - Q_pred).squaredNorm();
                b += (Q0 - Q_pred).squaredNorm();
            }
        }
    }
    if (a == 0 && b == 0) {
        return 1;
    }
    return 1 - sqrt(a / (b == 0 ? 1 : b));
}

void DistanceBarrierRBProblem::update_augmented_lagrangian(
    const Eigen::VectorXd& x)
{
    int ndof = PoseD::dim_to_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());
    int rot_ndof = PoseD::dim_to_rot_ndof(dim());

    double eta_q = compute_linear_augment_lagrangian_progress(x);
    double eta_Q = compute_angular_augment_lagrangian_progress(x);

    if (eta_q >= 0.999) {
        // Fix the kinematic DoF that have converged
        for (size_t i = 0; i < num_bodies(); i++) {
            if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
                is_dof_satisfied.segment(ndof * i, pos_ndof).setOnes();
            }
        }
    } else if (eta_q < 0.99 && linear_augmented_lagrangian_penalty < 1e8) {
        // Increase the κ_q
        linear_augmented_lagrangian_penalty *= 2;
    } else {
        // Increase the λ
        for (size_t i = 0, ki = 0; i < num_bodies(); i++) {
            if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
                linear_augmented_lagrangian_multiplier.segment(
                    ki * pos_ndof, pos_ndof) -=
                    linear_augmented_lagrangian_penalty
                    * sqrt(m_assembler[i].mass)
                    * (x.segment(ndof * i, pos_ndof)
                       - x_pred.segment(ndof * i, pos_ndof));
                ki++;
            }
        }
    }

    if (eta_Q >= 0.999) {
        // Fix the kinematic DoF that have converged
        for (size_t i = 0; i < num_bodies(); i++) {
            if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
                is_dof_satisfied.segment(ndof * i + pos_ndof, rot_ndof)
                    .setOnes();
            }
        }
    } else if (eta_Q < 0.99 && angular_augmented_lagrangian_penalty < 1e8) {
        // Increase the κ_Q
        angular_augmented_lagrangian_penalty *= 2;
    } else {
        // Increase the Λ
        for (size_t i = 0, ki = 0; i < num_bodies(); i++) {
            if (m_assembler[i].type == RigidBodyType::KINEMATIC) {
                size_t ri = ndof * i + pos_ndof;
                if (dim() == 2) {
                    angular_augmented_lagrangian_multiplier.middleRows(
                        ki * rot_ndof, rot_ndof) -=
                        angular_augmented_lagrangian_penalty
                        * sqrt(m_assembler[i].moment_of_inertia[0])
                        * (x.segment(ri, rot_ndof)
                           - x_pred.segment(ri, rot_ndof));
                } else {
                    auto Q_pred = construct_rotation_matrix(
                        VectorMax3d(x_pred.segment(ri, rot_ndof)));
                    auto Q = construct_rotation_matrix(
                        VectorMax3d(x.segment(ri, rot_ndof)));
                    angular_augmented_lagrangian_multiplier.middleRows(
                        rot_ndof * ki, rot_ndof) -=
                        angular_augmented_lagrangian_penalty * (Q - Q_pred)
                        * compute_Jsqrt(m_assembler[i].moment_of_inertia);
                }
                ki++;
            }
        }
    }

    if (eta_q < 0.999 || eta_Q < 0.999) {
        spdlog::info(
            "updated augmented Lagrangian "
            "κ_q={:g} κ_Q={:g} ||λ||∞={:g} ||Λ||∞={:g} η_q={:g} η_Q={:g}",
            linear_augmented_lagrangian_penalty,
            angular_augmented_lagrangian_penalty,
            linear_augmented_lagrangian_multiplier.lpNorm<Eigen::Infinity>(),
            angular_augmented_lagrangian_multiplier.lpNorm<Eigen::Infinity>(),
            eta_q, eta_Q);
    }
}

bool DistanceBarrierRBProblem::are_equality_constraints_satisfied(
    const Eigen::VectorXd& x) const
{
    if (m_assembler.count_kinematic_bodies()) {
        return compute_linear_augment_lagrangian_progress(x) >= 0.999
            && compute_angular_augment_lagrangian_progress(x) >= 0.999;
    }
    return true;
}

OptimizationResults DistanceBarrierRBProblem::solve_constraints()
{
    OptimizationResults opt_result;
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

        PosesD poses = this->dofs_to_poses(opt_result.x);

        Constraints collision_constraints;
        m_constraint.construct_constraint_set(
            m_assembler, poses, collision_constraints);
        update_friction_constraints(collision_constraints, poses);

        Eigen::VectorXd grad_Ex, grad_Bx, grad_Dx;
        compute_energy_term(opt_result.x, grad_Ex);
        compute_barrier_term(opt_result.x, collision_constraints, grad_Bx);
        compute_friction_term(opt_result.x, grad_Dx);

        Eigen::VectorXd tmp = grad_Ex + barrier_stiffness() * grad_Bx + grad_Dx;
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
            "Ending friction solve early because newton solve {:d} failed!", i);
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

    const double h = timestep();

    // update final pose
    // -------------------------------------
    m_assembler.set_rb_poses(this->dofs_to_poses(x));
    PosesD poses_q1 = m_assembler.rb_poses_t1();

    // Update the velocities
    // This need to be done AFTER updating poses
    for (RigidBody& rb : m_assembler.m_rbs) {
        if (rb.type != RigidBodyType::DYNAMIC) {
            continue;
        }

        // Linear update
        switch (body_energy_integration_method) {
        case IMPLICIT_EULER:
            rb.velocity.position =
                (rb.pose.position - rb.pose_prev.position) / h;
            break;

        case IMPLICIT_NEWMARK:
            rb.velocity.position =
                2 * (rb.pose.position - rb.pose_prev.position) / h
                - rb.velocity.position;
            rb.acceleration.position =
                2 * (rb.velocity.position - rb.velocity_prev.position) / h
                - rb.acceleration.position;
            break;
        case STABILIZED_NEWMARK: {
            rb.velocity.position =
                2 * (rb.pose.position - rb.pose_prev.position) / h
                - rb.velocity.position;
            // q̃ = q⁰ + hv⁰+ ¼h²(g + m⁻¹f + a⁰)
            VectorMax3d pos_tilde = rb.pose_prev.position
                + h
                    * (rb.velocity_prev.position
                       + h / 4.0
                           * (gravity + rb.force.position / rb.mass
                              + rb.acceleration.position));
            rb.acceleration.position =
                4 * (rb.pose.position - pos_tilde) / (h * h) + gravity
                + rb.force.position / rb.mass;
            break;
        }
        }

        // Angular update
        if (dim() == 2) {
            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                rb.velocity.rotation =
                    (rb.pose.rotation - rb.pose_prev.rotation) / h;
                break;

            case IMPLICIT_NEWMARK:
            case STABILIZED_NEWMARK:
                rb.velocity.rotation =
                    2 * (rb.pose.rotation - rb.pose_prev.rotation) / h
                    - rb.velocity.rotation;
                rb.acceleration.rotation =
                    2 * (rb.velocity.rotation - rb.velocity_prev.rotation) / h
                    - rb.acceleration.rotation;
                break;
            }
        } else {
            // Compute the rotation R s.t.
            // R * Rᵗ = Rᵗ⁺¹ → R = Rᵗ⁺¹(Rᵗ)ᵀ
            Eigen::Matrix3d R = rb.pose.construct_rotation_matrix()
                * rb.pose_prev.construct_rotation_matrix().transpose();
            // TODO: Make sure we did not loose momentum do to π modulus
            // ω = rotation_vector(R)
            Eigen::AngleAxisd omega(R);
            rb.velocity.rotation =
                omega.angle() / timestep() * rb.R0.transpose() * omega.axis();

            Eigen::Matrix3d Q = rb.pose.construct_rotation_matrix();
            Eigen::Matrix3d Q_prev = rb.pose_prev.construct_rotation_matrix();
            // Q̇ = Q[ω]
            // Q̇ᵗ = (Qᵗ - Qᵗ⁻¹) / h
            // Eigen::Matrix3d omega_hat = Q.transpose() * Qdot;
            // std::cout << omega_hat << std::endl << std::endl;
            // rb.velocity.rotation.x() = omega_hat(2, 1);
            // rb.velocity.rotation.y() = omega_hat(0, 2);
            // rb.velocity.rotation.z() = omega_hat(1, 0);
            // rb.velocity.rotation = omega.angle() / h
            //      * rb.R0.transpose() * omega.axis();

            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                rb.Qdot = (Q - Q_prev) / h;
                break;
            case IMPLICIT_NEWMARK: {
                auto Qdot_prev = rb.Qdot;
                rb.Qdot = 2 * (Q - Q_prev) / h - rb.Qdot;
                rb.Qddot = 2 * (rb.Qdot - Qdot_prev) / h - rb.Qddot;
                break;
            }
            case STABILIZED_NEWMARK: {
                auto Qdot_prev = rb.Qdot;
                rb.Qdot = 2 * (Q - Q_prev) / h - rb.Qdot;
                // auto Jinv = compute_Jinv(rb.moment_of_inertia);
                // Eigen::Matrix3d Tau =
                //     Q_prev.transpose() * Hat(rb.force.rotation);
                Eigen::Matrix3d Q_tilde = Q_prev
                    + h
                        * (Qdot_prev
                           + h / 4.0
                               * (
                                     // Tau * Jinv +
                                     rb.Qddot));
                rb.Qddot = 4 * (Q - Q_tilde) / (h * h); //+ Tau * Jinv;
                break;
            }
            }
        }

        rb.velocity.zero_dof(rb.is_dof_fixed, rb.R0);
        rb.acceleration.zero_dof(rb.is_dof_fixed, rb.R0);
    }

    if (do_intersection_check) {
        // Check for intersections instead of collision along the entire
        // step. We only guarentee a piecewise collision-free trajectory.
        // return detect_collisions(poses_t0, poses_q1,
        // CollisionCheck::EXACT);
        return detect_intersections(poses_q1);
    }
    return false;
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
    double Ex = compute_energy_term(x, grad, hess, compute_grad, compute_hess);
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
    Constraints constraints;
    m_constraint.construct_constraint_set(
        m_assembler, this->dofs_to_poses(x), constraints);

    spdlog::debug(
        "problem={} num_vertex_vertex_constraint={:d} "
        "num_edge_vertex_constraints={:d} num_edge_edge_constraints={:d} "
        "num_face_vertex_constraints={:d}",
        name(), constraints.vv_constraints.size(),
        constraints.ev_constraints.size(), constraints.ee_constraints.size(),
        constraints.fv_constraints.size());

    Eigen::VectorXd grad_Bx;
    Eigen::SparseMatrix<double> hess_Bx;
    double Bx = compute_barrier_term(
        x, constraints, grad_Bx, hess_Bx, compute_grad, compute_hess);

    // D(x) is the friction potential (Equation 15 in the IPC paper)
    Eigen::VectorXd grad_Dx;
    Eigen::SparseMatrix<double> hess_Dx;
    double Dx =
        compute_friction_term(x, grad_Dx, hess_Dx, compute_grad, compute_hess);

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

    int ndof = PoseD::dim_to_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());
    int rot_ndof = PoseD::dim_to_rot_ndof(dim());

    Eigen::VectorXd energies = Eigen::VectorXd::Zero(num_bodies());
    if (compute_grad) {
        grad.setZero(x.size());
    }
    tbb::concurrent_vector<Eigen::Triplet<double>> hess_triplets;
    if (compute_hess) {
        // Hessian is a block diagonal with (ndof x ndof) blocks
        hess_triplets.reserve(num_bodies() * ndof * ndof);
    }

    const std::vector<PoseD> poses = this->dofs_to_poses(x);
    assert(poses.size() == num_bodies());

    // tbb::parallel_for(size_t(0), poses.size(), [&](size_t i) {
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), poses.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            // Activate autodiff with the correct number of variables
            Diff::activate(ndof);

            for (long i = range.begin(); i != range.end(); ++i) {
                const PoseD& pose = poses[i];
                const RigidBody& body = m_assembler[i];

                // Do not compute the body energy for static and kinematic
                // bodies
                if (body.type != RigidBodyType::DYNAMIC) {
                    continue;
                }

                VectorMax6d gradi;

                if (compute_hess) {
                    // Initialize autodiff variables
                    Pose<Diff::DDouble2> pose_diff(Diff::d2vars(0, pose.dof()));

                    Diff::DDouble2 dExi = compute_body_energy<Diff::DDouble2>(
                        body, pose_diff,
                        grad_barrier_t0.segment(i * ndof, ndof));

                    energies[i] = dExi.getValue();
                    gradi = dExi.getGradient();
                    MatrixMax6d hessi = dExi.getHessian();

                    // Project dense block to make assembled matrix PSD
                    hessi.topLeftCorner(pos_ndof, pos_ndof) = project_to_psd(
                        MatrixMax3d(hessi.topLeftCorner(pos_ndof, pos_ndof)));
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
                    Pose<Diff::DDouble1> pose_diff(Diff::d1vars(0, pose.dof()));

                    Diff::DDouble1 dExi = compute_body_energy<Diff::DDouble1>(
                        body, pose_diff,
                        grad_barrier_t0.segment(i * ndof, ndof));

                    energies[i] = dExi.getValue();
                    gradi = dExi.getGradient();
                } else {
                    energies[i] = compute_body_energy<double>(
                        body, pose, grad_barrier_t0.segment(i * ndof, ndof));
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

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
    if (!is_checking_derivative) {
        is_checking_derivative = true;
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
        is_checking_derivative = false;
    }
#endif

    return energies.sum();
}

// Compute the energy term for a single rigid body
template <typename T>
T DistanceBarrierRBProblem::compute_body_energy(
    const RigidBody& body,
    const Pose<T>& pose,
    const VectorMax6d& grad_barrier_t0)
{
    // NOTE: t0 suffix indicates the current value not the inital value
    double h = timestep();

    T energy(0.0);

    // Linear energy
    if (!body.is_dof_fixed.head(pose.pos_ndof()).all()) {
        VectorMax3<T> q = pose.position;
        const VectorMax3d& q_t0 = body.pose.position;
        const VectorMax3d& qdot_t0 = body.velocity.position;
        VectorMax3d qddot_t0 = gravity + body.force.position / body.mass;
        switch (body_energy_integration_method) {
        case IMPLICIT_EULER:
            break;
        case IMPLICIT_NEWMARK:
        case STABILIZED_NEWMARK:
            qddot_t0 += body.acceleration.position;
            qddot_t0 *= 0.25;
            break;
        }

        // ½mqᵀq - mqᵀ(qᵗ + h(q̇ᵗ + h(g + f/m + ½∇B(qᵗ)/m)))
        energy += 0.5 * body.mass * q.dot(q)
            - body.mass * q.dot(q_t0 + h * (qdot_t0 + h * qddot_t0));
    }

    // Rotational energy
    if (!body.is_dof_fixed.tail(pose.rot_ndof()).all()) {
        if (dim() == 3) {
            Matrix3<T> Q = pose.construct_rotation_matrix();
            Eigen::Matrix3d Q_t0 = body.pose.construct_rotation_matrix();
            // Eigen::Matrix3d Qdot_t0 =
            //     Q_t0 * Hat(body.velocity.rotation);
            Eigen::Matrix3d Qdot_t0 = body.Qdot;

            DiagonalMatrix3d J = compute_J(body.moment_of_inertia);

            // Transform the world space torque into body space
            Eigen::Matrix3d Qddot_t0;
            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                Qddot_t0.setZero();
                break;
            case IMPLICIT_NEWMARK:
            case STABILIZED_NEWMARK:
                Qddot_t0 = 0.25 * body.Qddot;
                break;
            }

            // ½tr(QJQᵀ) - tr(Q(J(Qᵗ + hQ̇ᵗ + h²Aᵗ)ᵀ + h²[τ]))
            energy += 0.5 * (Q * J * Q.transpose()).trace();
            energy -=
                (Q * J * (Q_t0 + h * (Qdot_t0 + h * Qddot_t0)).transpose())
                    .trace();
            // Transform the world space torque into body space
            Eigen::Matrix3d Tau = Q_t0.transpose() * Hat(body.force.rotation);
            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                energy += h * h * (Q * Tau).trace();
                break;
            case IMPLICIT_NEWMARK:
            case STABILIZED_NEWMARK:
                energy += 0.25 * h * h * (Q * Tau).trace();
                break;
            }
        } else {
            assert(pose.rot_ndof() == 1);
            T theta = pose.rotation[0];
            double theta_t0 = body.pose.rotation[0];
            double theta_dot_t0 = body.velocity.rotation[0];
            // θ̈ = α + τ/I
            double I = body.moment_of_inertia[0];
            double theta_ddot_t0 = body.force.rotation[0] / I;
            switch (body_energy_integration_method) {
            case IMPLICIT_EULER:
                break;
            case IMPLICIT_NEWMARK:
            case STABILIZED_NEWMARK:
                theta_ddot_t0 += body.acceleration.rotation[0];
                theta_ddot_t0 *= 0.25;
                break;
            }

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
    int ndof = PoseD::dim_to_ndof(dim());
    int pos_ndof = PoseD::dim_to_pos_ndof(dim());
    int rot_ndof = PoseD::dim_to_rot_ndof(dim());
    size_t num_kinematic_bodies = m_assembler.count_kinematic_bodies();

    double potential = 0;
    if (compute_grad) {
        grad.setZero(x.size());
    }
    std::vector<Eigen::Triplet<double>> hess_triplets;
    if (compute_hess) {
        hess.resize(x.size(), x.size());
        hess_triplets.reserve(num_kinematic_bodies * ndof);
    }

    bool all_kinematic_dof_satisfied = true;
    for (size_t i = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type == RigidBodyType::KINEMATIC
            && !is_dof_satisfied.segment(ndof * i, ndof).all()) {
            all_kinematic_dof_satisfied = false;
            break;
        }
    }

    if (all_kinematic_dof_satisfied) {
        return potential;
    }

    PROFILE_POINT("DistanceBarrierRBProblem::compute_augmented_lagrangian");
    PROFILE_START();

    // Compute the linear AL potential
    const double& kappa_q = linear_augmented_lagrangian_penalty;
    for (size_t i = 0, ki = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type != RigidBodyType::KINEMATIC) {
            continue;
        }

        double m = m_assembler[i].mass;
        const auto& lambda = linear_augmented_lagrangian_multiplier.segment(
            ki * pos_ndof, pos_ndof);

        const auto& q = x.segment(i * ndof, pos_ndof);
        const auto& q_pred = x_pred.segment(i * ndof, pos_ndof);

        potential += kappa_q / 2 * m * (q - q_pred).squaredNorm()
            - sqrt(m) * lambda.dot(q - q_pred);
        if (compute_grad) {
            grad.segment(i * ndof, pos_ndof) =
                kappa_q * m * (q - q_pred) - sqrt(m) * lambda;
        }
        if (compute_hess) {
            for (int j = 0; j < pos_ndof; j++) {
                hess_triplets.emplace_back(
                    ndof * i + j, ndof * i + j, kappa_q * m);
            }
        }

        ki++;
    }

    typedef AutodiffType<Eigen::Dynamic, /*maxN=*/3> Diff;
    Diff::activate(rot_ndof);

    // Compute the angular AL potential
    const double& kappa_Q = angular_augmented_lagrangian_penalty;
    for (size_t i = 0, ki = 0; i < num_bodies(); i++) {
        if (m_assembler[i].type != RigidBodyType::KINEMATIC) {
            continue;
        }

        const VectorMax3d& moment_of_inertia = m_assembler[i].moment_of_inertia;
        MatrixMax3d lambda = angular_augmented_lagrangian_multiplier.middleRows(
            rot_ndof * ki, rot_ndof);

        VectorMax3d theta = x.segment(i * ndof + pos_ndof, rot_ndof);
        VectorMax3d theta_pred = x_pred.segment(i * ndof + pos_ndof, rot_ndof);

        if (dim() == 2) {
            double I = moment_of_inertia(0);
            double Isqrt = sqrt(I);

            potential += kappa_Q / 2 * I * (theta - theta_pred).squaredNorm()
                - (Isqrt * lambda.transpose() * (theta - theta_pred)).trace();

            if (compute_grad) {
                grad.segment(i * ndof + pos_ndof, rot_ndof) =
                    kappa_Q * I * (theta - theta_pred) - Isqrt * lambda;
            }
            if (compute_hess) {
                for (int j = 0; j < rot_ndof; j++) {
                    hess_triplets.emplace_back(
                        ndof * i + pos_ndof + j, ndof * i + pos_ndof + j,
                        kappa_Q * I);
                }
            }
        } else {
            VectorMax3<Diff::DDouble2> theta_diff = Diff::d2vars(0, theta);

            DiagonalMatrix3d J = compute_J(moment_of_inertia);
            DiagonalMatrix3d Jsqrt = compute_Jsqrt(moment_of_inertia);

            const auto& Q = construct_rotation_matrix(theta_diff);
            const auto& Q_pred = construct_rotation_matrix(theta_pred);

            Diff::DDouble2 dAL = kappa_Q / 2
                    * ((Q - Q_pred) * J * (Q - Q_pred).transpose()).trace()
                - (lambda.transpose() * (Q - Q_pred) * Jsqrt).trace();

            potential += dAL.getValue();
            if (compute_grad) {
                grad.segment(i * ndof + pos_ndof, rot_ndof) = dAL.getGradient();
            }
            if (compute_hess) {
                Eigen::Matrix3d H = dAL.getHessian();
                for (int hi = 0; hi < H.rows(); hi++) {
                    for (int hj = 0; hj < H.cols(); hj++) {
                        hess_triplets.emplace_back(
                            ndof * i + pos_ndof + hi, ndof * i + pos_ndof + hj,
                            H(hi, hj));
                    }
                }
            }
        }

        ki++;
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

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
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
    PosesD poses = this->dofs_to_poses(x);
    Constraints constraints;
    m_constraint.construct_constraint_set(m_assembler, poses, constraints);
    num_constraints = constraints.num_constraints();

    m_num_contacts = std::max(m_num_contacts, num_constraints);

    spdlog::debug(
        "problem={} num_vertex_vertex_constraint={:d} "
        "num_edge_vertex_constraints={:d} num_edge_edge_constraints={:d} "
        "num_face_vertex_constraints={:d}",
        name(), constraints.vv_constraints.size(),
        constraints.ev_constraints.size(), constraints.ee_constraints.size(),
        constraints.fv_constraints.size());

    double Bx = compute_barrier_term(
        x, constraints, grad, hess, compute_grad, compute_hess);

    return Bx;
}

// Convert from a local hessian to the triplets in the global hessian
template <typename DerivedLocalGradient>
void local_gradient_to_global(
    const Eigen::MatrixBase<DerivedLocalGradient>& local_gradient,
    const std::array<long, 2>& body_ids,
    int ndof,
    Eigen::VectorXd& grad)
{
    assert(local_gradient.size() == 2 * ndof);
    for (int b_i = 0; b_i < body_ids.size(); b_i++) {
        grad.segment(ndof * body_ids[b_i], ndof) +=
            local_gradient.segment(ndof * b_i, ndof);
    }
}

template <typename DerivedLocalHessian>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<DerivedLocalHessian>& local_hessian,
    const std::array<long, 2>& body_ids,
    int ndof,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    assert(local_hessian.rows() == 2 * ndof);
    assert(local_hessian.cols() == 2 * ndof);
    for (int b_i = 0; b_i < body_ids.size(); b_i++) {
        for (int b_j = 0; b_j < body_ids.size(); b_j++) {
            for (int dof_i = 0; dof_i < ndof; dof_i++) {
                for (int dof_j = 0; dof_j < ndof; dof_j++) {
                    double v =
                        local_hessian(ndof * b_i + dof_i, ndof * b_j + dof_j);
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
    const VectorMax12d& grad_f,
    const Eigen::MatrixXd& jac_V,
    const MatrixMax12d& hess_f,
    const Eigen::MatrixXd& hess_V,
    const std::vector<long>& vertex_ids,
    const std::vector<uint8_t>& local_body_ids,
    const std::array<long, 2>& body_ids,
    const int dim,
    Eigen::VectorXd& grad,
    std::vector<Eigen::Triplet<double>>& hess_triplets,
    bool compute_grad,
    bool compute_hess)
{
    if (!compute_grad && !compute_hess) {
        return;
    }

    // PROFILE_POINT("apply_chain_rule");
    // PROFILE_START();

    const int rb_ndof = PoseD::dim_to_ndof(dim);

    if (compute_grad) {
        // jac_Vi ∈ R^{4n × 2m}
        VectorMax12d local_grad = VectorMax12d::Zero(2 * rb_ndof);
        for (int i = 0; i < vertex_ids.size(); i++) {
            local_grad.segment(rb_ndof * local_body_ids[i], rb_ndof) +=
                jac_V.middleRows(vertex_ids[i] * dim, dim).transpose()
                * grad_f.segment(i * dim, dim);
        }

        local_gradient_to_global(local_grad, body_ids, rb_ndof, grad);
    }

    if (compute_hess) {
        // jac_Vi ∈ R^{4n × 2m}
        MatrixMax12d jac_Vi =
            MatrixMax12d::Zero(vertex_ids.size() * dim, 2 * rb_ndof);
        for (int i = 0; i < vertex_ids.size(); i++) {
            jac_Vi.block(i * dim, local_body_ids[i] * rb_ndof, dim, rb_ndof) =
                jac_V.middleRows(vertex_ids[i] * dim, dim);
        }

        // hess ∈ R^{2m × 2m}
        MatrixMax12d hess = jac_Vi.transpose() * hess_f * jac_Vi;
        for (int i = 0; i < vertex_ids.size(); i++) {
            for (int j = 0; j < dim; j++) {
                // Off diagaonal blocks are all zero because the derivative
                // of a vertex of body A with body B is zero.
                hess.block(
                    local_body_ids[i] * rb_ndof, local_body_ids[i] * rb_ndof,
                    rb_ndof, rb_ndof) +=
                    hess_V.middleRows(
                        rb_ndof * (vertex_ids[i] * dim + j), rb_ndof)
                    * grad_f[i * dim + j];
            }
        }

        hess = project_to_psd(hess);

        local_hessian_to_global_triplets(
            hess, body_ids, rb_ndof, hess_triplets);
    }

    // PROFILE_END();
}

struct PotentialStorage {
    PotentialStorage() {}
    PotentialStorage(size_t nvars) { gradient.setZero(nvars); }
    double potential = 0;
    Eigen::VectorXd gradient;
    std::vector<Eigen::Triplet<double>> hessian_triplets;
};
typedef tbb::enumerable_thread_specific<PotentialStorage>
    ThreadSpecificPotentials;

double merge_derivative_storage(
    const ThreadSpecificPotentials& potentials,
    size_t nvars,
    Eigen::VectorXd& grad,
    Eigen::SparseMatrix<double>& hess,
    bool compute_grad,
    bool compute_hess)
{
    PROFILE_POINT("merge_derivative_storage");
    PROFILE_START();

    if (compute_grad) {
        grad.setZero(nvars);
    }
    if (compute_hess) {
        hess.resize(nvars, nvars);
    }

    double potential = 0;
    for (const auto& p : potentials) {
        potential += p.potential;

        if (compute_grad) {
            grad += p.gradient;
        }

        if (compute_hess) {
            Eigen::SparseMatrix<double> p_hess(nvars, nvars);
            p_hess.setFromTriplets(
                p.hessian_triplets.begin(), p.hessian_triplets.end());
            hess += p_hess;
        }
    }

    PROFILE_END();

    return potential;
}

double DistanceBarrierRBProblem::compute_barrier_term(
    const Eigen::VectorXd& x,
    const Constraints& constraints,
    Eigen::VectorXd& grad,
    Eigen::SparseMatrix<double>& hess,
    bool compute_grad,
    bool compute_hess)
{
    if (constraints.size() == 0) {
        grad.setZero(x.size());
        hess.resize(x.size(), x.size());
        return 0;
    }

    PROFILE_POINT("DistanceBarrierRBProblem::compute_barrier_term");
    // WARNING: PROFILE_POINTs are not thread safe
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierRBProblem::compute_barrier_term:value",
    //     COMPUTE_BARRIER_VAL);
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierRBProblem::compute_barrier_term:gradient",
    //     COMPUTE_BARRIER_GRAD);
    // NAMED_PROFILE_POINT(
    //     "DistanceBarrierRBProblem::compute_barrier_term:hessian",
    //     COMPUTE_BARRIER_HESS);

    PROFILE_START();

    int rb_ndof = PoseD::dim_to_ndof(dim());

    // Compute V(x)
    Eigen::MatrixXd jac_V, hess_V;
    Eigen::MatrixXd V = m_assembler.world_vertices_diff(
        x, jac_V, hess_V, compute_grad || compute_hess, compute_hess);

    double dhat = barrier_activation_distance();

    ThreadSpecificPotentials thread_storage(x.size());
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), constraints.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            // Get references to the local derivative storage
            auto& local_storage = thread_storage.local();
            auto& potential = local_storage.potential;
            auto& local_grad = local_storage.gradient;
            auto& hess_triplets = local_storage.hessian_triplets;

            for (size_t ci = range.begin(); ci != range.end(); ++ci) {
                const auto& constraint = constraints[ci];

                // PROFILE_START(COMPUTE_BARRIER_VAL);
                potential +=
                    constraint.compute_potential(V, edges(), faces(), dhat);
                // PROFILE_START(COMPUTE_BARRIER_VAL);

                VectorMax12d grad_B;
                if (compute_grad || compute_hess) {
                    // PROFILE_START(COMPUTE_BARRIER_GRAD);
                    grad_B = constraint.compute_potential_gradient(
                        V, edges(), faces(), dhat);
                    // PROFILE_END(COMPUTE_BARRIER_GRAD);
                }

                MatrixMax12d hess_B;
                if (compute_hess) {
                    // PROFILE_START(COMPUTE_BARRIER_HESS);
                    hess_B = constraint.compute_potential_hessian(
                        V, edges(), faces(), dhat,
                        /*project_hessian_to_psd=*/false);
                    // PROFILE_END(COMPUTE_BARRIER_HESS);
                }

                apply_chain_rule(
                    grad_B, jac_V, hess_B, hess_V,
                    constraint.vertex_indices(edges(), faces()),
                    vertex_local_body_ids(constraints, ci),
                    body_ids(m_assembler, constraints, ci), dim(), local_grad,
                    hess_triplets, compute_grad, compute_hess);
            }
        });

    double potential = merge_derivative_storage(
        thread_storage, x.size(), grad, hess, compute_grad, compute_hess);

    PROFILE_END();

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
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

// NAMED_PROFILE_POINT(
//     "DistanceBarrierRBProblem::compute_friction_potential:value",
//     COMPUTE_FRICTION_VAL);
// NAMED_PROFILE_POINT(
//     "DistanceBarrierRBProblem::compute_friction_potential:gradient",
//     COMPUTE_FRICTION_GRAD);
// NAMED_PROFILE_POINT(
//     "DistanceBarrierRBProblem::compute_friction_potential:hessian",
//     COMPUTE_FRICTION_HESS);
template <typename RigidBodyConstraint, typename FrictionConstraint>
double DistanceBarrierRBProblem::compute_friction_potential(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXd& jac_V,
    const Eigen::MatrixXd& hess_V,
    const FrictionConstraint& constraint,
    Eigen::VectorXd& grad,
    std::vector<Eigen::Triplet<double>>& hess_triplets,
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

    int rb_ndof = PoseD::dim_to_ndof(dim());

    double epsv_times_h = static_friction_speed_bound * timestep();

    // PROFILE_START(COMPUTE_FRICTION_VAL);
    double Dx = constraint.compute_potential(U, edges(), faces(), epsv_times_h);
    // PROFILE_END(COMPUTE_FRICTION_VAL);

    VectorMax12d grad_D;
    if (compute_grad || compute_hess) {
        // PROFILE_START(COMPUTE_FRICTION_GRAD);
        grad_D = constraint.compute_potential_gradient(
            U, edges(), faces(), epsv_times_h);
        // PROFILE_END(COMPUTE_FRICTION_GRAD);
    }

    MatrixMax12d hess_D;
    if (compute_hess) {
        // PROFILE_START(COMPUTE_FRICTION_HESS);
        hess_D = constraint.compute_potential_hessian(
            U, edges(), faces(), epsv_times_h,
            /*project_hessian_to_psd=*/false);
        // PROFILE_END(COMPUTE_FRICTION_HESS);
    }

    RigidBodyConstraint rbc(m_assembler, constraint);
    apply_chain_rule(
        grad_D, jac_V, hess_D, hess_V,
        constraint.vertex_indices(edges(), faces()),
        rbc.vertex_local_body_ids(), rbc.body_ids(), dim(), //
        grad, hess_triplets, compute_grad, compute_hess);

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
        hess.resize(x.size(), x.size());
        return 0;
    }

    PROFILE_POINT("DistanceBarrierRBProblem::compute_friction_term");
    PROFILE_START();

    int rb_ndof = PoseD::dim_to_ndof(dim());

    // Compute V(x)
    Eigen::MatrixXd jac_V, hess_V;
    Eigen::MatrixXd V1 = m_assembler.world_vertices_diff(
        x, jac_V, hess_V, compute_grad || compute_hess, compute_hess);

    NAMED_PROFILE_POINT(
        "DistanceBarrierRBProblem::compute_friction_term:displacement",
        DISPLACEMENT);
    PROFILE_START(DISPLACEMENT);
    // absolute linear dislacement of each point
    Eigen::MatrixXd U = V1 - m_assembler.world_vertices(poses_t0);
    PROFILE_END(DISPLACEMENT);

    ThreadSpecificPotentials thread_storage(x.size());
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraints.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            // Get references to the local derivative storage
            auto& local_storage = thread_storage.local();
            auto& potential = local_storage.potential;
            auto& local_grad = local_storage.gradient;
            auto& hess_triplets = local_storage.hessian_triplets;

            for (size_t ci = range.begin(); ci != range.end(); ++ci) {
                size_t local_ci = ci;

                if (local_ci < friction_constraints.vv_constraints.size()) {
                    potential += compute_friction_potential<
                        RigidBodyVertexVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.vv_constraints[local_ci],
                        local_grad, hess_triplets, compute_grad, compute_hess);
                    continue;
                }

                local_ci -= friction_constraints.vv_constraints.size();
                if (local_ci < friction_constraints.ev_constraints.size()) {
                    potential += compute_friction_potential<
                        RigidBodyEdgeVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.ev_constraints[local_ci],
                        local_grad, hess_triplets, compute_grad, compute_hess);
                    continue;
                }

                local_ci -= friction_constraints.ev_constraints.size();
                if (local_ci < friction_constraints.ee_constraints.size()) {
                    potential +=
                        compute_friction_potential<RigidBodyEdgeEdgeConstraint>(
                            U, jac_V, hess_V,
                            friction_constraints.ee_constraints[local_ci],
                            local_grad, hess_triplets, compute_grad,
                            compute_hess);
                    continue;
                }

                local_ci -= friction_constraints.ee_constraints.size();
                assert(local_ci < friction_constraints.fv_constraints.size());
                potential +=
                    compute_friction_potential<RigidBodyFaceVertexConstraint>(
                        U, jac_V, hess_V,
                        friction_constraints.fv_constraints[local_ci],
                        local_grad, hess_triplets, compute_grad, compute_hess);
            }
        });

    double potential = merge_derivative_storage(
        thread_storage, x.size(), grad, hess, compute_grad, compute_hess);

    PROFILE_END();

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
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

    return potential;
}

///////////////////////////////////////////////////////////////////////////

double DistanceBarrierRBProblem::compute_min_distance() const
{
    double min_distance = m_constraint.compute_minimum_distance(
        m_assembler, m_assembler.rb_poses());
    return std::isfinite(min_distance) ? min_distance : -1;
}

double
DistanceBarrierRBProblem::compute_min_distance(const Eigen::VectorXd& x) const
{
    PosesD poses = this->dofs_to_poses(x);
    double min_distance =
        m_constraint.compute_minimum_distance(m_assembler, poses);
    return std::isfinite(min_distance) ? min_distance : -1;
}

bool DistanceBarrierRBProblem::has_collisions(
    const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j)
{
    PosesD poses_i = this->dofs_to_poses(x_i);
    PosesD poses_j = this->dofs_to_poses(x_j);
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

    PosesD poses_i = this->dofs_to_poses(x_i);
    PosesD poses_j = this->dofs_to_poses(x_j);
    double earliest_toi =
        m_constraint.compute_earliest_toi(m_assembler, poses_i, poses_j);
    m_had_collisions |= earliest_toi <= 1;
    return earliest_toi;
}

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
// The following functions are used exclusivly to check that the
// gradient and hessian match a finite difference version.

void DistanceBarrierRBProblem::check_barrier_gradient(
    const Eigen::VectorXd& x,
    const Constraints& constraints,
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
    const Constraints& constraints,
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
    // hess_approx = project_to_psd(hess_approx);
    if (!fd::compare_jacobian(hess, hess_approx, Constants::FINITE_DIFF_TEST)) {
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
    auto f = [&](const Eigen::VectorXd& x) { return compute_friction_term(x); };
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
    hess_approx = project_to_psd(hess_approx);
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

    Eigen::MatrixXd hess_autodiff = project_to_psd(f_diff.getHessian());

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
        spdlog::error("finite gradient check failed for augmented lagrangian");
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
        spdlog::error("finite hessian check failed for augmented lagrangian");
    }
}
#endif

} // namespace ipc::rigid
