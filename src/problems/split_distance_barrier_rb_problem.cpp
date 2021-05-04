#include "split_distance_barrier_rb_problem.hpp"

#include <tbb/parallel_sort.h>

#include <finitediff.hpp>

#include <constants.hpp>
#include <logger.hpp>
#include <profiler.hpp>
#include <solvers/solver_factory.hpp>
#include <time_stepper/time_stepper_factory.hpp>

namespace ipc::rigid {

SplitDistanceBarrierRBProblem::SplitDistanceBarrierRBProblem()
    : DistanceBarrierRBProblem()
{
}

bool SplitDistanceBarrierRBProblem::settings(const nlohmann::json& params)
{
    auto time_stepper_name =
        params["rigid_body_problem"]["time_stepper"].get<std::string>();
    m_time_stepper = time_stepper_name == "default"
        ? TimeStepperFactory::factory().get_default_time_stepper(dim())
        : TimeStepperFactory::factory().get_time_stepper(time_stepper_name);
    return DistanceBarrierRBProblem::settings(params);
}

nlohmann::json SplitDistanceBarrierRBProblem::settings() const
{
    nlohmann::json json = DistanceBarrierRBProblem::settings();
    json["time_stepper"] = m_time_stepper->name();
    return json;
}

////////////////////////////////////////////////////////////
// Rigid Body Problem

void SplitDistanceBarrierRBProblem::simulation_step(
    bool& had_collision, bool& _has_intersections, bool solve_collision)
{
    // Take an unconstrained time-step
    m_time_stepper->step(m_assembler, gravity, timestep());

    update_dof();

    had_collision =
        detect_collisions(poses_t0, poses_t1, CollisionCheck::CONSERVATIVE);

    // Check if minimum distance is violated
    min_distance = compute_min_distance(this->poses_to_dofs(poses_t1));
    if (min_distance >= 0) {
        spdlog::info("candidate_step min_distance={:.8e}", min_distance);

        // our constraint is really d > min_d, we want to run the
        // optimization when we end the step with small distances
        if (min_distance <= barrier_activation_distance()) {
            had_collision = true;
        }
    }

    if (solve_collision && had_collision) {
        update_constraints();
        OptimizationResults result = solve_constraints();
        _has_intersections = take_step(result.x);
    } else {
        _has_intersections = detect_intersections(this->poses_t1);
    }
}

void SplitDistanceBarrierRBProblem::update_dof()
{
    DistanceBarrierRBProblem::update_dof();
    poses_t1 = m_assembler.rb_poses_t1();
}

void SplitDistanceBarrierRBProblem::update_constraints()
{
    DistanceBarrierRBProblem::update_constraints();

    // Compute the collisions
    original_impacts.clear();

    // TODO: original_impacts is only used in the 2D restitution code
    // if this changes remove this if
    if (dim() == 2) {
        constraint().construct_collision_set(
            m_assembler, poses_t0, poses_t1, original_impacts);

        tbb::parallel_sort(
            original_impacts.ev_impacts.begin(),
            original_impacts.ev_impacts.end(),
            compare_impacts_by_time<EdgeVertexImpact>);
        tbb::parallel_sort(
            original_impacts.ee_impacts.begin(),
            original_impacts.ee_impacts.end(),
            compare_impacts_by_time<EdgeEdgeImpact>);
        tbb::parallel_sort(
            original_impacts.fv_impacts.begin(),
            original_impacts.fv_impacts.end(),
            compare_impacts_by_time<FaceVertexImpact>);
    }
}

bool SplitDistanceBarrierRBProblem::take_step(const Eigen::VectorXd& dof)
{
    if (coefficient_restitution < 0) {
        return DistanceBarrierRBProblem::take_step(dof);
    }

    // This need to be done BEFORE updating poses
    // -------------------------------------
    solve_velocities();

    // update final pose
    // -------------------------------------
    m_assembler.set_rb_poses(this->dofs_to_poses(dof));
    PosesD poses_q1 = m_assembler.rb_poses_t1();

    // Check for intersections instead of collision along the entire
    // step. We only guarentee a piecewise collision-free trajectory.
    // return detect_collisions(poses_t0, poses_q1, CollisionCheck::EXACT);
    return detect_intersections(poses_q1);
}

void SplitDistanceBarrierRBProblem::solve_velocities()
{
    if (dim() != 2) {
        throw NotImplementedError(
            "RigidBodyProblem::solve_velocities() not implmented for 3D!");
    }

    // precompute normal directions (since velocities will change i
    // can't do it after
    Eigen::MatrixXd normals(original_impacts.ev_impacts.size(), 2);

    for (long i = 0; i < normals.rows(); ++i) {
        const EdgeVertexImpact& ev_impact =
            original_impacts.ev_impacts[size_t(i)];

        double toi = ev_impact.time;

        const long edge_id = ev_impact.edge_index;
        const long a_id = ev_impact.vertex_index;
        const int b0_id = m_assembler.m_edges.coeff(edge_id, 0);
        const int b1_id = m_assembler.m_edges.coeff(edge_id, 1);

        const size_t body_B_id =
            size_t(m_assembler.m_vertex_to_body_map(b0_id));
        bool is_oriented = m_assembler[body_B_id].is_oriented;

        Eigen::Vector2d n_toi;
        VectorMax3d e_toi; // edge vector at toi
        switch (constraint().trajectory_type) {
        case TrajectoryType::LINEAR: {
            // Use linearized trajectories
            VectorMax3d e_v0_t0 = m_assembler.world_vertex(poses_t0, b0_id);
            VectorMax3d e_v0_t1 = m_assembler.world_vertex(poses_t1, b0_id);
            VectorMax3d e_v0_toi = (e_v0_t1 - e_v0_t0) * toi + e_v0_t0;

            VectorMax3d e_v1_t0 = m_assembler.world_vertex(poses_t0, b1_id);
            VectorMax3d e_v1_t1 = m_assembler.world_vertex(poses_t1, b1_id);
            VectorMax3d e_v1_toi = (e_v1_t1 - e_v1_t0) * toi + e_v1_t0;

            e_toi = e_v1_toi - e_v0_toi;
        } break;

        case TrajectoryType::PIECEWISE_LINEAR:
        case TrajectoryType::RIGID:
        case TrajectoryType::REDON: {
            // Use nonlinear trajectory
            long edge_body_id = m_assembler.edge_id_to_body_id(edge_id);

            PoseD pose_toi = PoseD::interpolate(
                poses_t0[edge_body_id], poses_t1[edge_body_id], toi);

            e_toi = m_assembler.world_vertex(pose_toi, b1_id)
                - m_assembler.world_vertex(pose_toi, b0_id);
        } break;
        }
        n_toi << -e_toi(1), e_toi(0); // 90deg ccw rotation
        n_toi.normalize();

        if (is_oriented) {
            n_toi = -n_toi;
        } else {
            // check normal points towards A
            Eigen::Vector2d va = m_assembler.world_vertex(poses_t0, a_id);
            Eigen::Vector2d vb = m_assembler.world_vertex(poses_t0, b0_id);
            if ((va - vb).transpose() * n_toi <= 0.0) {
                n_toi *= -1;
            }
        }
        normals.row(i) = n_toi.transpose();
    }

#ifndef NDEBUG
    double prev_toi = -1;
#endif
    for (long i = 0; i < normals.rows(); ++i) {
        const EdgeVertexImpact& ev_impact =
            original_impacts.ev_impacts[size_t(i)];

        double toi = ev_impact.time;
#ifndef NDEBUG
        assert(prev_toi <= toi);
        prev_toi = toi;
#endif
        double alpha = ev_impact.alpha;

        // global ids of the vertices
        const long edge_id = ev_impact.edge_index;
        const long a_id = ev_impact.vertex_index;
        const int b0_id = m_assembler.m_edges.coeff(edge_id, 0);
        const int b1_id = m_assembler.m_edges.coeff(edge_id, 1);

        const size_t body_A_id = size_t(m_assembler.m_vertex_to_body_map(a_id));
        const size_t body_B_id =
            size_t(m_assembler.m_vertex_to_body_map(b0_id));

        // local (rigid body) ids of the vertices
        const long r_A_id = a_id - m_assembler.m_body_vertex_id[body_A_id];
        const long r_B0_id = b0_id - m_assembler.m_body_vertex_id[body_B_id];
        const long r_B1_id = b1_id - m_assembler.m_body_vertex_id[body_B_id];

        auto& body_A = m_assembler[body_A_id];
        auto& body_B = m_assembler[body_B_id];

        // The velocities of the center of mass at the time of
        // collision!!
        PoseD vel_A_prev =
            PoseD::interpolate(body_A.velocity_prev, body_A.velocity, toi);
        PoseD vel_B_prev =
            PoseD::interpolate(body_B.velocity_prev, body_B.velocity, toi);

        // The masss
        const double inv_m_A =
            body_A.is_dof_fixed.head(dim()).any() ? 0.0 : 1.0 / body_A.mass;
        const double inv_m_B =
            body_B.is_dof_fixed.head(dim()).any() ? 0.0 : 1.0 / body_B.mass;

        if (dim() != 2) {
            throw NotImplementedError("Resitution not implemented in 3D yet!");
        }

        // The moment of inertia
        int rot_ndof = PoseD::dim_to_rot_ndof(dim());
        const Eigen::MatrixXd inv_I_A = body_A.is_dof_fixed.tail(rot_ndof).any()
            ? Eigen::MatrixXd::Zero(rot_ndof, rot_ndof)
            : Eigen::MatrixXd(
                  body_A.moment_of_inertia.cwiseInverse().asDiagonal());
        const Eigen::MatrixXd inv_I_B = body_B.is_dof_fixed.tail(rot_ndof).any()
            ? Eigen::MatrixXd::Zero(rot_ndof, rot_ndof)
            : Eigen::MatrixXd(
                  body_B.moment_of_inertia.cwiseInverse().asDiagonal());

        // Vectors from the center of mass to the collision point
        // (90deg rotation counter clockwise)
        //
        // (1) first get vertices position wrt rigid bodies
        const Eigen::Vector2d r0_A = body_A.vertices.row(r_A_id);
        const Eigen::Vector2d r0_B0 =
            body_B.vertices.row(r_B0_id); // edge vertex 0
        const Eigen::Vector2d r0_B1 =
            body_B.vertices.row(r_B1_id); // edge vertex 1
        const Eigen::Vector2d r0_B = r0_B0 + alpha * (r0_B1 - r0_B0);

        // (2) and the angular displacement at time of collision
        const PoseD pose_Atoi =
            PoseD::interpolate(body_A.pose_prev, body_A.pose, toi);
        const PoseD pose_Btoi =
            PoseD::interpolate(body_B.pose_prev, body_B.pose, toi);
        // (3) then the vectors are given by r = R(\theta_{t})*r_0
        Eigen::Matrix2d one_hat = Hat(1.0);
        const Eigen::VectorXd r_Aperp_toi =
            pose_Atoi.construct_rotation_matrix() * one_hat * r0_A;
        const Eigen::VectorXd r_Bperp_toi =
            pose_Btoi.construct_rotation_matrix() * one_hat * r0_B;

        // The collision point velocities BEFORE collision
        const Eigen::VectorXd v_Aprev =
            vel_A_prev.position + vel_A_prev.rotation(0) * r_Aperp_toi;
        const Eigen::VectorXd v_Bprev =
            vel_B_prev.position + vel_B_prev.rotation(0) * r_Bperp_toi;

        // The relative veolicity magnitud BEFORE collision
        const Eigen::VectorXd& n_toi = normals.row(i);

        const double vrel_prev_toi = (v_Aprev - v_Bprev).transpose() * n_toi;
        if (vrel_prev_toi >= 0.0) {
            continue;
        }
        // solve for the impulses
        const double nr_A_toi = n_toi.transpose() * r_Aperp_toi;
        const double nr_B_toi = n_toi.transpose() * r_Bperp_toi;
        const double K =
            (inv_m_A + inv_m_B //
             + inv_I_A(0) * nr_A_toi * nr_A_toi
             + inv_I_B(0) * nr_B_toi * nr_B_toi);

        const double j =
            -1.0 / K * (1.0 + coefficient_restitution) * vrel_prev_toi;

        // update
        Eigen::Vector2d V_A_delta = inv_m_A * j * n_toi;
        Eigen::Vector2d V_B_delta = -inv_m_B * j * n_toi;
        double w_A_delta = inv_I_A(0) * j * nr_A_toi;
        double w_B_delta = -inv_I_B(0) * j * nr_B_toi;

        if (!(body_A.is_dof_fixed[0] || body_A.is_dof_fixed[1])) {
            body_A.velocity.position = vel_A_prev.position + V_A_delta;
        }

        if (!(body_B.is_dof_fixed[0] || body_B.is_dof_fixed[1])) {
            body_B.velocity.position = vel_B_prev.position + V_B_delta;
        }

        if (!body_A.is_dof_fixed[2]) {
            body_A.velocity.rotation(0) = vel_A_prev.rotation(0) + w_A_delta;
        }

        if (!body_B.is_dof_fixed[2]) {
            body_B.velocity.rotation(0) = vel_B_prev.rotation(0) + w_B_delta;
        }
    }
}

////////////////////////////////////////////////////////////
// Barrier Problem

// Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
double SplitDistanceBarrierRBProblem::compute_energy_term(
    const Eigen::VectorXd& x,
    Eigen::VectorXd& grad,
    Eigen::SparseMatrix<double>& hess,
    bool compute_grad,
    bool compute_hess)
{
    Eigen::VectorXd diff = x - this->poses_to_dofs(poses_t1);
    DiagonalMatrixXd M = m_assembler.m_rb_mass_matrix;

    if (compute_grad) {
        grad = M * diff;
#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
        Eigen::VectorXd grad_approx = eval_grad_energy_approx(*this, x);
        if (!fd::compare_gradient(grad, grad_approx)) {
            spdlog::error("finite gradient check failed for E(x)");
        }
#endif
    }

    if (compute_hess) {
        hess = SparseDiagonal(M.diagonal());
#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
        Eigen::MatrixXd hess_approx = eval_hess_energy_approx(*this, x);
        if (!fd::compare_jacobian(hess, hess_approx)) {
            spdlog::error("finite hessian check failed for E(x)");
        }
#endif
    }

    return 0.5 * diff.transpose() * M * diff;
}

} // namespace ipc::rigid
