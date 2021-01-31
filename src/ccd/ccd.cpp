#include "ccd.hpp"

#include <mutex>

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_invoke.h>

#include <ipc/ccd/ccd.hpp>
#include <ipc/friction/closest_point.hpp>

#include <ccd/linear/broad_phase.hpp>
#include <ccd/linear/edge_vertex_ccd.hpp>
#include <ccd/piecewise_linear/time_of_impact.hpp>
#include <ccd/redon/time_of_impact.hpp>
#include <ccd/rigid/broad_phase.hpp>
#include <ccd/rigid/time_of_impact.hpp>

// #define SAVE_CCD_QUERIES
#ifdef SAVE_CCD_QUERIES
#include <ccd/save_queries.hpp>
#endif

#include <profiler.hpp>

namespace ipc::rigid {

void detect_collisions(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Impacts& impacts,
    DetectionMethod method,
    TrajectoryType trajectory)
{
    assert(bodies.num_bodies() == poses_t0.size());
    assert(poses_t0.size() == poses_t1.size());

    // Do the broad phase by detecting candidate impacts
    Candidates candidates;
    detect_collision_candidates(
        bodies, poses_t0, poses_t1, collision_types, candidates, method,
        trajectory);

    // Do the narrow phase by detecting actual impacts from the candidate set
    detect_collisions_from_candidates(
        bodies, poses_t0, poses_t1, candidates, impacts, trajectory);
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collision_candidates(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const int collision_types,
    Candidates& candidates,
    DetectionMethod method,
    TrajectoryType trajectory,
    const double inflation_radius)
{
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    switch (trajectory) {
    case TrajectoryType::LINEAR:
        detect_collision_candidates_linear(
            bodies, poses_t0, poses_t1, collision_types, candidates, method,
            inflation_radius);
        break;
    case TrajectoryType::PIECEWISE_LINEAR:
    case TrajectoryType::RIGID:
    case TrajectoryType::REDON:
        detect_collision_candidates_rigid(
            bodies, poses_t0, poses_t1, collision_types, candidates, method,
            inflation_radius);
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Broad-Phase CCD
///////////////////////////////////////////////////////////////////////////////

void detect_collisions_from_candidates(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const Candidates& candidates,
    Impacts& impacts,
    TrajectoryType trajectory)
{
    PROFILE_POINT("collisions_detection__narrow_phase");
    PROFILE_START();

    std::mutex ev_impacts_mutex, ee_impacts_mutex, fv_impacts_mutex;

    auto ev_impact = [&](const EdgeVertexCandidate& ev_candidate) {
        double toi;
        bool is_colliding = edge_vertex_ccd(
            bodies, poses_t0, poses_t1, ev_candidate, toi, trajectory);
        if (is_colliding) {
            double alpha = edge_vertex_closest_point(
                bodies, poses_t0, poses_t1, ev_candidate, toi, trajectory);
            std::scoped_lock lock(ev_impacts_mutex);
            impacts.ev_impacts.emplace_back(
                toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index);
        }
    };

    auto ee_impact = [&](const EdgeEdgeCandidate& ee_candidate) {
        double toi;
        bool is_colliding = edge_edge_ccd(
            bodies, poses_t0, poses_t1, ee_candidate, toi, trajectory);
        if (is_colliding) {
            double alpha, beta;
            edge_edge_closest_point(
                bodies, poses_t0, poses_t1, ee_candidate, toi, alpha, beta,
                trajectory);
            std::scoped_lock lock(ee_impacts_mutex);
            impacts.ee_impacts.emplace_back(
                toi, ee_candidate.edge0_index, alpha, ee_candidate.edge1_index,
                beta);
        }
    };

    auto fv_impact = [&](const FaceVertexCandidate& fv_candidate) {
        double toi;
        bool is_colliding = face_vertex_ccd(
            bodies, poses_t0, poses_t1, fv_candidate, toi, trajectory);
        if (is_colliding) {
            double u, v;
            face_vertex_closest_point(
                bodies, poses_t0, poses_t1, fv_candidate, toi, u, v,
                trajectory);
            std::scoped_lock lock(fv_impacts_mutex);
            impacts.fv_impacts.emplace_back(
                toi, fv_candidate.face_index, u, v, fv_candidate.vertex_index);
        }
    };

    impacts.clear();
    tbb::parallel_invoke(
        [&] { tbb::parallel_for_each(candidates.ev_candidates, ev_impact); },
        [&] { tbb::parallel_for_each(candidates.ee_candidates, ee_impact); },
        [&] { tbb::parallel_for_each(candidates.fv_candidates, fv_impact); });

    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
bool edge_vertex_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi,
    double minimum_separation_distance)
{
    assert(bodies.dim() == 2);

#ifdef SAVE_CCD_QUERIES
    save_ccd_candidate(bodies, poses_t0, poses_t1, candidate);
#endif

    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long e0_id = bodyB.edges(edge_id, 0);
    long e1_id = bodyB.edges(edge_id, 1);

    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector2d v_t0 = bodyA.world_vertex(poseA_t0, vertex_id);
        Eigen::Vector2d v_t1 = bodyA.world_vertex(poseA_t1, vertex_id);

        Eigen::Vector2d e0_t0 = bodyB.world_vertex(poseB_t0, e0_id);
        Eigen::Vector2d e0_t1 = bodyB.world_vertex(poseB_t1, e0_id);

        Eigen::Vector2d e1_t0 = bodyB.world_vertex(poseB_t0, e1_id);
        Eigen::Vector2d e1_t1 = bodyB.world_vertex(poseB_t1, e1_id);

        // Expects the arguments as position and displacements
        return compute_edge_vertex_time_of_impact(
            e0_t0, e1_t0, v_t0, (e0_t1 - e0_t0).eval(), (e1_t1 - e1_t0).eval(),
            (v_t1 - v_t0).eval(), toi);
    }

    case TrajectoryType::PIECEWISE_LINEAR:
        return compute_piecewise_linear_edge_vertex_time_of_impact(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            edge_id, toi, earliest_toi, minimum_separation_distance);

    case TrajectoryType::RIGID:
        return compute_edge_vertex_time_of_impact(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            edge_id, toi, earliest_toi);

    case TrajectoryType::REDON:
        return compute_edge_vertex_time_of_impact_redon(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            edge_id, toi, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

bool edge_edge_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi,
    double minimum_separation_distance)
{
#ifdef SAVE_CCD_QUERIES
    save_ccd_candidate(bodies, poses_t0, poses_t1, candidate);
#endif

    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long ea0_id = bodyA.edges(edgeA_id, 0);
    long ea1_id = bodyA.edges(edgeA_id, 1);
    long eb0_id = bodyB.edges(edgeB_id, 0);
    long eb1_id = bodyB.edges(edgeB_id, 1);

    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector3d ea0_t0 = bodyA.world_vertex(poseA_t0, ea0_id);
        Eigen::Vector3d ea0_t1 = bodyA.world_vertex(poseA_t1, ea0_id);

        Eigen::Vector3d ea1_t0 = bodyA.world_vertex(poseA_t0, ea1_id);
        Eigen::Vector3d ea1_t1 = bodyA.world_vertex(poseA_t1, ea1_id);

        Eigen::Vector3d eb0_t0 = bodyB.world_vertex(poseB_t0, eb0_id);
        Eigen::Vector3d eb0_t1 = bodyB.world_vertex(poseB_t1, eb0_id);

        Eigen::Vector3d eb1_t0 = bodyB.world_vertex(poseB_t0, eb1_id);
        Eigen::Vector3d eb1_t1 = bodyB.world_vertex(poseB_t1, eb1_id);

        // TODO: Check if edges are parallel
        bool result = ipc::edge_edge_ccd(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi,
            /*conservative_rescaling=*/0.8);
        assert(!result || toi > 0);
        return result;
    }

    case TrajectoryType::PIECEWISE_LINEAR:
        return compute_piecewise_linear_edge_edge_time_of_impact(
            bodyA, poseA_t0, poseA_t1, edgeA_id, bodyB, poseB_t0, poseB_t1,
            edgeB_id, toi, earliest_toi, minimum_separation_distance);

    case TrajectoryType::RIGID:
        return compute_edge_edge_time_of_impact(
            bodyA, poseA_t0, poseA_t1, edgeA_id, bodyB, poseB_t0, poseB_t1,
            edgeB_id, toi, earliest_toi);

    case TrajectoryType::REDON:
        return compute_edge_edge_time_of_impact_redon(
            bodyA, poseA_t0, poseA_t1, edgeA_id, bodyB, poseB_t0, poseB_t1,
            edgeB_id, toi, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

bool face_vertex_ccd(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& candidate,
    double& toi,
    TrajectoryType trajectory,
    double earliest_toi,
    double minimum_separation_distance)
{
#ifdef SAVE_CCD_QUERIES
    save_ccd_candidate(bodies, poses_t0, poses_t1, candidate);
#endif

    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long f0_id = bodyB.faces(face_id, 0);
    long f1_id = bodyB.faces(face_id, 1);
    long f2_id = bodyB.faces(face_id, 2);

    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector3d v_t0 = bodyA.world_vertex(poseA_t0, vertex_id);
        Eigen::Vector3d v_t1 = bodyA.world_vertex(poseA_t1, vertex_id);

        Eigen::Vector3d f0_t0 = bodyB.world_vertex(poseB_t0, f0_id);
        Eigen::Vector3d f0_t1 = bodyB.world_vertex(poseB_t1, f0_id);

        Eigen::Vector3d f1_t0 = bodyB.world_vertex(poseB_t0, f1_id);
        Eigen::Vector3d f1_t1 = bodyB.world_vertex(poseB_t1, f1_id);

        Eigen::Vector3d f2_t0 = bodyB.world_vertex(poseB_t0, f2_id);
        Eigen::Vector3d f2_t1 = bodyB.world_vertex(poseB_t1, f2_id);

        bool result = ipc::point_triangle_ccd(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, toi,
            /*conservative_rescaling=*/0.8);
        assert(!result || toi > 0);
        return result;
    }

    case TrajectoryType::PIECEWISE_LINEAR:
        return compute_piecewise_linear_face_vertex_time_of_impact(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            face_id, toi, earliest_toi, minimum_separation_distance);

    case TrajectoryType::RIGID:
        return compute_face_vertex_time_of_impact(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            face_id, toi, earliest_toi);

    case TrajectoryType::REDON:
        return compute_face_vertex_time_of_impact_redon(
            bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
            face_id, toi, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

double edge_vertex_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& candidate,
    double toi,
    TrajectoryType trajectory)
{
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long e0_id = bodyB.edges(edge_id, 0);
    long e1_id = bodyB.edges(edge_id, 1);

    Eigen::Vector2d v, e0, e1;
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector2d v_t0 = bodyA.world_vertex(poseA_t0, vertex_id);
        Eigen::Vector2d v_t1 = bodyA.world_vertex(poseA_t1, vertex_id);
        v = (v_t1 - v_t0) * toi + v_t0;

        Eigen::Vector2d e0_t0 = bodyB.world_vertex(poseB_t0, e0_id);
        Eigen::Vector2d e0_t1 = bodyB.world_vertex(poseB_t1, e0_id);
        e0 = (e0_t1 - e0_t0) * toi + e0_t0;

        Eigen::Vector2d e1_t0 = bodyB.world_vertex(poseB_t0, e1_id);
        Eigen::Vector2d e1_t1 = bodyB.world_vertex(poseB_t1, e1_id);
        e1 = (e1_t1 - e1_t0) * toi + e1_t0;

        break;
    }

    case TrajectoryType::PIECEWISE_LINEAR:
    case TrajectoryType::RIGID:
    case TrajectoryType::REDON: {
        // Compute the poses at time toi
        PoseD poseA_toi = PoseD::interpolate(poseA_t0, poseA_t1, toi);
        PoseD poseB_toi = PoseD::interpolate(poseB_t0, poseB_t1, toi);
        v = bodyA.world_vertex(poseA_toi, vertex_id);
        e0 = bodyB.world_vertex(poseB_toi, e0_id);
        e1 = bodyB.world_vertex(poseB_toi, e1_id);
        break;
    }
    }
    return ipc::point_edge_closest_point(v, e0, e1);
}

void edge_edge_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& candidate,
    double toi,
    double& alpha,
    double& beta,
    TrajectoryType trajectory)
{
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long ea0_id = bodyA.edges(edgeA_id, 0);
    long ea1_id = bodyA.edges(edgeA_id, 1);
    long eb0_id = bodyB.edges(edgeB_id, 0);
    long eb1_id = bodyB.edges(edgeB_id, 1);

    Eigen::Vector3d ea0, ea1, eb0, eb1;
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector3d ea0_t0 = bodyA.world_vertex(poseA_t0, ea0_id);
        Eigen::Vector3d ea0_t1 = bodyA.world_vertex(poseA_t1, ea0_id);
        ea0 = (ea0_t1 - ea0_t0) * toi + ea0_t0;

        Eigen::Vector3d ea1_t0 = bodyA.world_vertex(poseA_t0, ea1_id);
        Eigen::Vector3d ea1_t1 = bodyA.world_vertex(poseA_t1, ea1_id);
        ea1 = (ea1_t1 - ea1_t0) * toi + ea1_t0;

        Eigen::Vector3d eb0_t0 = bodyB.world_vertex(poseB_t0, eb0_id);
        Eigen::Vector3d eb0_t1 = bodyB.world_vertex(poseB_t1, eb0_id);
        eb0 = (eb0_t1 - eb0_t0) * toi + eb0_t0;

        Eigen::Vector3d eb1_t0 = bodyB.world_vertex(poseB_t0, eb1_id);
        Eigen::Vector3d eb1_t1 = bodyB.world_vertex(poseB_t1, eb1_id);
        eb1 = (eb1_t1 - eb1_t0) * toi + eb1_t0;

        break;
    }

    case TrajectoryType::PIECEWISE_LINEAR:
    case TrajectoryType::RIGID:
    case TrajectoryType::REDON: {
        // Compute the poses at time toi
        PoseD poseA_toi = PoseD::interpolate(poseA_t0, poseA_t1, toi);
        PoseD poseB_toi = PoseD::interpolate(poseB_t0, poseB_t1, toi);

        ea0 = bodyA.world_vertex(poseA_toi, bodyA.edges(edgeA_id, 0));
        ea1 = bodyA.world_vertex(poseA_toi, bodyA.edges(edgeA_id, 1));

        eb0 = bodyB.world_vertex(poseB_toi, bodyB.edges(edgeB_id, 0));
        eb1 = bodyB.world_vertex(poseB_toi, bodyB.edges(edgeB_id, 1));
        break;
    }
    }
    Eigen::Vector2d alpha_beta =
        ipc::edge_edge_closest_point(ea0, ea1, eb0, eb1);
    alpha = alpha_beta[0];
    beta = alpha_beta[1];
}

void face_vertex_closest_point(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& candidate,
    double toi,
    double& u,
    double& v,
    TrajectoryType trajectory)
{
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);
    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];
    const PoseD& poseA_t0 = poses_t0[bodyA_id];
    const PoseD& poseA_t1 = poses_t1[bodyA_id];
    const PoseD& poseB_t0 = poses_t0[bodyB_id];
    const PoseD& poseB_t1 = poses_t1[bodyB_id];
    long f0_id = bodyB.faces(face_id, 0);
    long f1_id = bodyB.faces(face_id, 1);
    long f2_id = bodyB.faces(face_id, 2);

    Eigen::Vector3d p, f0, f1, f2;
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::Vector3d p_t0 = bodyA.world_vertex(poseA_t0, vertex_id);
        Eigen::Vector3d p_t1 = bodyA.world_vertex(poseA_t1, vertex_id);
        p = (p_t1 - p_t0) * toi + p_t0;

        Eigen::Vector3d f0_t0 = bodyB.world_vertex(poseB_t0, f0_id);
        Eigen::Vector3d f0_t1 = bodyB.world_vertex(poseB_t1, f0_id);
        f0 = (f0_t1 - f0_t0) * toi + f0_t0;

        Eigen::Vector3d f1_t0 = bodyB.world_vertex(poseB_t0, f1_id);
        Eigen::Vector3d f1_t1 = bodyB.world_vertex(poseB_t1, f1_id);
        f1 = (f1_t1 - f1_t0) * toi + f1_t0;

        Eigen::Vector3d f2_t0 = bodyB.world_vertex(poseB_t0, f2_id);
        Eigen::Vector3d f2_t1 = bodyB.world_vertex(poseB_t1, f2_id);
        f2 = (f2_t1 - f2_t0) * toi + f2_t0;
        break;
    }

    case TrajectoryType::PIECEWISE_LINEAR:
    case TrajectoryType::RIGID:
    case TrajectoryType::REDON: {
        // Compute the poses at time toi
        PoseD poseA_toi = PoseD::interpolate(poseA_t0, poseA_t1, toi);
        PoseD poseB_toi = PoseD::interpolate(poseB_t0, poseB_t1, toi);

        // Get the world vertex of the point at time t
        p = bodyA.world_vertex(poseA_toi, vertex_id);
        // Get the world vertex of the face at time t
        f0 = bodyB.world_vertex(poseB_toi, f0_id);
        f1 = bodyB.world_vertex(poseB_toi, f1_id);
        f2 = bodyB.world_vertex(poseB_toi, f2_id);
        break;
    }
    }
    Eigen::Vector2d uv = ipc::point_triangle_closest_point(p, f0, f1, f2);
    u = uv[0];
    v = uv[1];
}

} // namespace ipc::rigid
