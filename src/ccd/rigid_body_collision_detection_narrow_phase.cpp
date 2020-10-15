#include "rigid_body_collision_detection.hpp"

#include <cmath>

#include <igl/barycentric_coordinates.h>

#include <ipc/friction/closest_point.hpp>

#include <ccd/rigid_body_time_of_impact.hpp>
#include <constants.hpp>
#include <geometry/intersection.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const Candidates& candidates,
    ConcurrentImpacts& impacts,
    TrajectoryType trajectory)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses_t0);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses_t1);
        detect_collisions_from_candidates(
            V_t0, V_t1, bodies.m_edges, bodies.m_faces, candidates, impacts);
        return;
    }
    assert(trajectory == TrajectoryType::SCREWING);

    NAMED_PROFILE_POINT(
        "rigid_collisions_detection__narrow_phase", NARROW_PHASE);
    PROFILE_START(NARROW_PHASE);

    auto detect_ev_collision = [&](const EdgeVertexCandidate& ev_candidate) {
        double toi, alpha;
        bool is_colliding = detect_edge_vertex_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, ev_candidate, toi, alpha, trajectory);
        if (is_colliding) {
            impacts.ev_impacts.emplace_back(
                toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index);
        }
    };

    auto detect_ee_collision = [&](const EdgeEdgeCandidate& ee_candidate) {
        double toi, edge0_alpha, edge1_alpha;
        bool is_colliding = detect_edge_edge_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, ee_candidate, toi, edge0_alpha,
            edge1_alpha, trajectory);
        if (is_colliding) {
            impacts.ee_impacts.emplace_back(
                toi, ee_candidate.edge0_index, edge0_alpha,
                ee_candidate.edge1_index, edge1_alpha);
        }
    };

    auto detect_fv_collision = [&](const FaceVertexCandidate& fv_candidate) {
        double toi, u, v;
        bool is_colliding = detect_face_vertex_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, fv_candidate, toi, u, v, trajectory);
        if (is_colliding) {
            impacts.fv_impacts.emplace_back(
                toi, fv_candidate.face_index, u, v, fv_candidate.vertex_index);
        }
    };

    tbb::parallel_invoke(
        [&] {
            impacts.ev_impacts.clear();
            tbb::parallel_for_each(
                candidates.ev_candidates, detect_ev_collision);
        },

        [&] {
            impacts.ee_impacts.clear();
            tbb::parallel_for_each(
                candidates.ee_candidates, detect_ee_collision);
        },

        [&] {
            impacts.fv_impacts.clear();
            tbb::parallel_for_each(
                candidates.fv_candidates, detect_fv_collision);
        });

    PROFILE_END(NARROW_PHASE);
}

// Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const EdgeVertexCandidate& candidate,
    double& toi,
    double& alpha,
    TrajectoryType trajectory,
    double earliest_toi)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::VectorX3d edge_vertex0_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge_index, 0));
        Eigen::VectorX3d edge_vertex0_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge_index, 0));

        Eigen::VectorX3d edge_vertex1_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge_index, 1));
        Eigen::VectorX3d edge_vertex1_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge_index, 1));

        Eigen::VectorX3d vertex_t0 =
            bodies.world_vertex(poses_t0, candidate.vertex_index);
        Eigen::VectorX3d vertex_t1 =
            bodies.world_vertex(poses_t1, candidate.vertex_index);

        return detect_edge_vertex_collisions_narrow_phase(
            edge_vertex0_t0, edge_vertex1_t0, vertex_t0, //
            edge_vertex0_t1, edge_vertex1_t1, vertex_t1, toi, alpha);
    }

    // TODO: Make this thread safe
    // PROFILE_POINT("detect_edge_vertex_collisions_narrow_phase");
    // PROFILE_START();

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    bool is_colliding = compute_edge_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses_t0[bodyA_id], poses_t1[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses_t0[bodyB_id],
        poses_t1[bodyB_id], edge_id, toi, earliest_toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA = physics::Pose<double>::interpolate(
            poses_t0[bodyA_id], poses_t1[bodyA_id], toi);
        physics::Pose<double> poseB = physics::Pose<double>::interpolate(
            poses_t0[bodyB_id], poses_t1[bodyB_id], toi);

        // Get the world vertex of the point at time t
        Eigen::VectorX3d vertex =
            bodies.m_rbs[bodyA_id].world_vertex(poseA, vertex_id);
        // Get the world vertex of the edge at time t
        Eigen::VectorX3d edge_vertex0 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 0));
        Eigen::VectorX3d edge_vertex1 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 1));

        // Compute the impact parameter along the edge
        alpha =
            ipc::point_edge_closest_point(vertex, edge_vertex0, edge_vertex1);
        // assert(
        //     alpha > -Constants::PARAMETER_ASSERTION_TOL
        //     && alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
    }

    // PROFILE_END();

    return is_colliding;
}

bool detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const EdgeEdgeCandidate& candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    TrajectoryType trajectory,
    double earliest_toi)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::VectorX3d edge0_vertex0_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge0_index, 0));
        Eigen::VectorX3d edge0_vertex0_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge0_index, 0));

        Eigen::VectorX3d edge0_vertex1_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge0_index, 1));
        Eigen::VectorX3d edge0_vertex1_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge0_index, 1));

        Eigen::VectorX3d edge1_vertex0_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge1_index, 0));
        Eigen::VectorX3d edge1_vertex0_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge1_index, 0));

        Eigen::VectorX3d edge1_vertex1_t0 = bodies.world_vertex(
            poses_t0, bodies.m_edges(candidate.edge1_index, 1));
        Eigen::VectorX3d edge1_vertex1_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge1_index, 1));

        return detect_edge_edge_collisions_narrow_phase(
            edge0_vertex0_t0, edge0_vertex1_t0, edge1_vertex0_t0,
            edge1_vertex1_t0, edge0_vertex0_t1, edge0_vertex1_t1,
            edge1_vertex0_t1, edge1_vertex1_t1, toi, edge0_alpha, edge1_alpha);
    }

    // TODO: Make this thread safe
    // PROFILE_POINT("detect_edge_edge_collisions_narrow_phase");
    // PROFILE_START();

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);

    bool is_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses_t0[bodyA_id], poses_t1[bodyA_id],
        edgeA_id, bodies.m_rbs[bodyB_id], poses_t0[bodyB_id],
        poses_t1[bodyB_id], edgeB_id, toi, earliest_toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA = physics::Pose<double>::interpolate(
            poses_t0[bodyA_id], poses_t1[bodyA_id], toi);
        physics::Pose<double> poseB = physics::Pose<double>::interpolate(
            poses_t0[bodyB_id], poses_t1[bodyB_id], toi);

        Eigen::Vector3d edgeA_vertex0_toi = bodies.m_rbs[bodyA_id].world_vertex(
            poseA, bodies.m_rbs[bodyA_id].edges(edgeA_id, 0));
        Eigen::Vector3d edgeA_vertex1_toi = bodies.m_rbs[bodyA_id].world_vertex(
            poseA, bodies.m_rbs[bodyA_id].edges(edgeA_id, 1));

        Eigen::Vector3d edgeB_vertex0_toi = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edgeB_id, 0));
        Eigen::Vector3d edgeB_vertex1_toi = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edgeB_id, 1));

        Eigen::Vector2d alphas = ipc::edge_edge_closest_point(
            edgeA_vertex0_toi, edgeA_vertex1_toi, //
            edgeB_vertex0_toi, edgeB_vertex1_toi);
        edge0_alpha = alphas(0);
        edge1_alpha = alphas(1);
        // assert(
        //     edge0_alpha > -Constants::PARAMETER_ASSERTION_TOL
        //     && edge0_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
        // assert(
        //     edge1_alpha > -Constants::PARAMETER_ASSERTION_TOL
        //     && edge1_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
    }

    // PROFILE_END();

    return is_colliding;
}

bool detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const FaceVertexCandidate& candidate,
    double& toi,
    double& u,
    double& v,
    TrajectoryType trajectory,
    double earliest_toi)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::VectorX3d face_vertex0_t0 = bodies.world_vertex(
            poses_t0, bodies.m_faces(candidate.face_index, 0));
        Eigen::VectorX3d face_vertex0_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 0));

        Eigen::VectorX3d face_vertex1_t0 = bodies.world_vertex(
            poses_t0, bodies.m_faces(candidate.face_index, 1));
        Eigen::VectorX3d face_vertex1_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 1));

        Eigen::VectorX3d face_vertex2_t0 = bodies.world_vertex(
            poses_t0, bodies.m_faces(candidate.face_index, 2));
        Eigen::VectorX3d face_vertex2_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 2));

        Eigen::VectorX3d vertex_t0 =
            bodies.world_vertex(poses_t0, candidate.vertex_index);
        Eigen::VectorX3d vertex_t1 =
            bodies.world_vertex(poses_t1, candidate.vertex_index);

        return detect_face_vertex_collisions_narrow_phase(
            face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, vertex_t0,
            face_vertex0_t1, face_vertex1_t1, face_vertex2_t1, vertex_t1, //
            toi, u, v);
    }

    // TODO: Make this thread safe
    // PROFILE_POINT("detect_face_vertex_collisions_narrow_phase");
    // PROFILE_START();

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);

    bool is_colliding = compute_face_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses_t0[bodyA_id], poses_t1[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses_t0[bodyB_id],
        poses_t1[bodyB_id], face_id, toi, earliest_toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA = physics::Pose<double>::interpolate(
            poses_t0[bodyA_id], poses_t1[bodyA_id], toi);
        physics::Pose<double> poseB = physics::Pose<double>::interpolate(
            poses_t0[bodyB_id], poses_t1[bodyB_id], toi);

        // Get the world vertex of the point at time t
        Eigen::RowVector3d vertex_toi =
            bodies.m_rbs[bodyA_id].world_vertex(poseA, vertex_id);
        // Get the world vertex of the face at time t
        Eigen::RowVector3d triangle_vertex0_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 0));
        Eigen::RowVector3d triangle_vertex1_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 1));
        Eigen::RowVector3d triangle_vertex2_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 2));

        // TODO: Consider moving this computation to an as needed basis
        Eigen::RowVector3d coords;
        igl::barycentric_coordinates(
            vertex_toi, triangle_vertex0_toi, triangle_vertex1_toi,
            triangle_vertex2_toi, coords);
        u = coords(0);
        v = coords(1);
        double w = coords(2);
        assert(
            u + v + w < 1 + Constants::PARAMETER_ASSERTION_TOL
            && u + v + w > 1 - Constants::PARAMETER_ASSERTION_TOL);
    }

    // PROFILE_END();

    return is_colliding;
}

} // namespace ccd
