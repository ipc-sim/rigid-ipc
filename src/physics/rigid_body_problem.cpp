#include "rigid_body_problem.hpp"

#include <iostream>

#include <tbb/parallel_for_each.h>

#include <finitediff.hpp>
#include <igl/PI.h>
#include <igl/predicates/segment_segment_intersect.h>

#include <ipc/utils/intersection.hpp>

#include <ccd/rigid/broad_phase.hpp>
#include <ccd/rigid/rigid_body_hash_grid.hpp>
#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <profiler.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

RigidBodyProblem::RigidBodyProblem()
    : coefficient_restitution(0)
    , coefficient_friction(0)
    , collision_eps(2)
    , m_timestep(0.01)
    , do_intersection_check(false)
{
    gravity.setZero(3);
}

bool RigidBodyProblem::settings(const nlohmann::json& params)
{
    collision_eps = params["collision_eps"];
    coefficient_restitution = params["coefficient_restitution"];
    coefficient_friction = params["coefficient_friction"];
    if (coefficient_friction < 0 || coefficient_friction > 1) {
        spdlog::warn(
            "Coefficient of friction (μ={:g}) is outside the standard "
            "[0, 1]",
            coefficient_friction);
    } else if (coefficient_friction == 0) {
        spdlog::info(
            "Disabling friction because coefficient of friction is zero");
    }

    std::vector<RigidBody> rbs;
    bool success = read_rb_scene(params, rbs);
    if (!success) {
        spdlog::error("Unable to read rigid body scene!");
        return false;
    }

    init(rbs);

    from_json(params["gravity"], gravity);
    assert(gravity.size() >= dim());
    gravity.conservativeResize(dim());

    do_intersection_check = params["do_intersection_check"];
    return true;
}

nlohmann::json RigidBodyProblem::settings() const
{
    nlohmann::json json;

    json["collision_eps"] = collision_eps;
    json["coefficient_restitution"] = coefficient_restitution;
    json["coefficient_friction"] = coefficient_friction;
    json["gravity"] = to_json(gravity);
    json["do_intersection_check"] = do_intersection_check;
    return json;
}

void RigidBodyProblem::init(const std::vector<RigidBody>& rbs)
{
    m_assembler.init(rbs);

    update_constraints();

    for (size_t i = 0; i < num_bodies(); ++i) {
        auto& rb = m_assembler[i];
        spdlog::info(
            "rb={:d} group_id={:d} mass={:g} inertia={}", i, rb.group_id,
            rb.mass, fmt_eigen(rb.moment_of_inertia));
    }
    spdlog::info("average_mass={:g}", m_assembler.average_mass);

    // Compute world diagonal
    Eigen::MatrixXd V = vertices();
    init_bbox_diagonal =
        (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
    spdlog::info("init_bbox_diagonal={:g}", init_bbox_diagonal);

    // Ensure the dimension of gravity matches the dimension of the problem.
    gravity = gravity.head(dim());

    if (detect_intersections(m_assembler.rb_poses_t1())) {
        spdlog::error("The initial state contains intersections!");
    } else {
        spdlog::info("no intersections found in initial state");
    }
}

nlohmann::json RigidBodyProblem::state() const
{
    nlohmann::json json;
    std::vector<nlohmann::json> rbs;
    Eigen::VectorXd p =
        Eigen::VectorXd::Zero(PoseD::dim_to_pos_ndof(dim())); // Linear momentum
    Eigen::VectorXd L = Eigen::VectorXd::Zero(
        PoseD::dim_to_rot_ndof(dim())); // Angular momentum
    double T = 0.0;                     // Kinetic energy
    double G = 0.0;                     // Potential energy

    for (auto& rb : m_assembler.m_rbs) {
        nlohmann::json jrb;
        jrb["position"] = to_json(Eigen::VectorXd(rb.pose.position));
        jrb["rotation"] = to_json(Eigen::VectorXd(rb.pose.rotation));
        jrb["linear_velocity"] = to_json(Eigen::VectorXd(rb.velocity.position));
        jrb["angular_velocity"] =
            to_json(Eigen::VectorXd(rb.velocity.rotation));
        if (dim() == 3) {
            jrb["Qdot"] = to_json(rb.Qdot);
            jrb["Qddot"] = to_json(rb.Qddot);
        }
        rbs.push_back(jrb);

        // momentum
        p += rb.mass * rb.velocity.position;
        L += rb.moment_of_inertia.asDiagonal() * rb.velocity.rotation;

        T += 0.5 * rb.mass * rb.velocity.position.squaredNorm();
        T += 0.5 * rb.velocity.rotation.transpose()
            * rb.moment_of_inertia.asDiagonal() * rb.velocity.rotation;

        if (!rb.is_dof_fixed[0] && !rb.is_dof_fixed[1]) {
            G -= rb.mass * gravity.dot(rb.pose.position);
        }
    }

    json["rigid_bodies"] = rbs;
    json["linear_momentum"] = to_json(p);
    json["angular_momentum"] = to_json(L);
    json["kinetic_energy"] = T;
    json["potential_energy"] = G;
    return json;
}

void RigidBodyProblem::state(const nlohmann::json& args)
{
    nlohmann::json json;
    auto& rbs = args["rigid_bodies"];
    assert(rbs.size() == num_bodies());
    size_t i = 0;
    for (auto& jrb : args["rigid_bodies"]) {
        from_json(jrb["position"], m_assembler[i].pose.position);
        from_json(jrb["rotation"], m_assembler[i].pose.rotation);
        from_json(jrb["linear_velocity"], m_assembler[i].velocity.position);
        from_json(jrb["angular_velocity"], m_assembler[i].velocity.rotation);
        if (dim() == 3) {
            if (jrb.contains("Qdot")) {
                from_json(jrb["Qdot"], m_assembler[i].Qdot);
            } else {
                spdlog::warn("Missing field \"Qdot\" in rigid body state!");
                m_assembler[i].Qdot.setZero();
            }
            if (jrb.contains("Qddot")) {
                from_json(jrb["Qddot"], m_assembler[i].Qddot);
            } else {
                spdlog::warn("Missing field \"Qddot\" in rigid body state!");
                m_assembler[i].Qddot.setZero();
            }
        }
        i++;
    }
}

void RigidBodyProblem::update_dof()
{
    poses_t0 = m_assembler.rb_poses_t0();
    x0 = this->poses_to_dofs(poses_t0);
    num_vars_ = x0.size();
}

void RigidBodyProblem::update_constraints()
{
    update_dof();
    constraint().initialize();
}

void RigidBodyProblem::init_solve()
{
    return solver().init_solve(starting_point());
}

OptimizationResults RigidBodyProblem::solve_constraints()
{
    return solver().solve(starting_point());
}

OptimizationResults RigidBodyProblem::step_solve()
{
    return solver().step_solve();
}

bool RigidBodyProblem::take_step(const Eigen::VectorXd& dof)
{
    ////////////////////////////////////////////////////////////////////////
    // WARNING: This only assumes an implicit euler velocity update. For
    // more updates look at the overridden version in
    // distance_barrier_rb_problem.
    ////////////////////////////////////////////////////////////////////////

    // update final pose
    // -------------------------------------
    m_assembler.set_rb_poses(this->dofs_to_poses(dof));
    PosesD poses_q1 = m_assembler.rb_poses_t1();

    // Update the velocities
    // This need to be done AFTER updating poses
    for (RigidBody& rb : m_assembler.m_rbs) {
        if (rb.type != RigidBodyType::DYNAMIC) {
            continue;
        }

        // Assume linear velocity through the time-step.
        rb.velocity.position =
            (rb.pose.position - rb.pose_prev.position) / timestep();
        if (dim() == 2) {
            rb.velocity.rotation =
                (rb.pose.rotation - rb.pose_prev.rotation) / timestep();
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

            // Q̇ = Q[ω]
            // Q̇ᵗ = (Qᵗ - Qᵗ⁻¹) / h
            Eigen::Matrix3d Q = rb.pose.construct_rotation_matrix();
            Eigen::Matrix3d Qdot =
                (Q - rb.pose_prev.construct_rotation_matrix()) / timestep();
            // Eigen::Matrix3d omega_hat = Q.transpose() * Qdot;
            // std::cout << omega_hat << std::endl << std::endl;
            // rb.velocity.rotation.x() = omega_hat(2, 1);
            // rb.velocity.rotation.y() = omega_hat(0, 2);
            // rb.velocity.rotation.z() = omega_hat(1, 0);
            rb.Qdot = Qdot;
        }
        rb.velocity.zero_dof(rb.is_dof_fixed, rb.R0);
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

bool RigidBodyProblem::detect_collisions(
    const PosesD& poses_q0,
    const PosesD& poses_q1,
    const CollisionCheck check_type) const
{
    Impacts impacts;

    double scale =
        check_type == CollisionCheck::EXACT ? 1.0 : (1.0 + collision_eps);
    PosesD scaled_pose_q1 = interpolate(poses_q0, poses_q1, scale);

    constraint().construct_collision_set(
        m_assembler, poses_q0, scaled_pose_q1, impacts);

    return impacts.size();
}

// Check if the geometry is intersecting
bool RigidBodyProblem::detect_intersections(const PosesD& poses) const
{
    if (num_bodies() <= 1) {
        return false;
    }

    PROFILE_POINT("RigidBodyProblem::detect_intersections");
    PROFILE_START();

    const Eigen::MatrixXd vertices = m_assembler.world_vertices(poses);
    const Eigen::MatrixXi& edges = this->edges();
    const Eigen::MatrixXi& faces = this->faces();

    bool is_intersecting = false;
    if (dim() == 2) { // Need to check segment-segment intersections in 2D
        assert(vertices.cols() == 2);

        double inflation_radius = 1e-8; // Conservative broad phase
        std::vector<std::pair<int, int>> close_bodies =
            m_assembler.close_bodies(poses, poses, inflation_radius);
        if (close_bodies.size() == 0) {
            PROFILE_END();
            return false;
        }

        RigidBodyHashGrid hashgrid;
        hashgrid.resize(m_assembler, poses, close_bodies, inflation_radius);
        hashgrid.addBodies(m_assembler, poses, close_bodies, inflation_radius);

        const Eigen::VectorXi& vertex_group_ids = group_ids();
        auto can_vertices_collide = [&vertex_group_ids](size_t vi, size_t vj) {
            return vertex_group_ids[vi] != vertex_group_ids[vj];
        };

        std::vector<EdgeEdgeCandidate> ee_candidates;
        hashgrid.getEdgeEdgePairs(edges, ee_candidates, can_vertices_collide);

        for (const EdgeEdgeCandidate& ee_candidate : ee_candidates) {
            if (igl::predicates::segment_segment_intersect(
                    vertices.row(edges(ee_candidate.edge0_index, 0)).head<2>(),
                    vertices.row(edges(ee_candidate.edge0_index, 1)).head<2>(),
                    vertices.row(edges(ee_candidate.edge1_index, 0)).head<2>(),
                    vertices.row(edges(ee_candidate.edge1_index, 1))
                        .head<2>())) {
                is_intersecting = true;
                break;
            }
        }
    } else { // Need to check segment-triangle intersections in 3D
        assert(dim() == 3);

        std::vector<EdgeFaceCandidate> ef_candidates;
        detect_intersection_candidates_rigid_bvh(
            m_assembler, poses, ef_candidates);

        for (const EdgeFaceCandidate& ef_candidate : ef_candidates) {
            if (is_edge_intersecting_triangle(
                    vertices.row(edges(ef_candidate.edge_index, 0)),
                    vertices.row(edges(ef_candidate.edge_index, 1)),
                    vertices.row(faces(ef_candidate.face_index, 0)),
                    vertices.row(faces(ef_candidate.face_index, 1)),
                    vertices.row(faces(ef_candidate.face_index, 2)))) {
                is_intersecting = true;
                break;
            }
        }
    }

    PROFILE_END();

    return is_intersecting;
}

} // namespace ipc::rigid
