#include "rigid_body_collision_detection.hpp"

#include <cmath>

#include <ccd/rigid_body_time_of_impact.hpp>
#include <constants.hpp>
#include <geometry/barycentric_coordinates.hpp>
#include <geometry/intersection.hpp>
#include <geometry/projection.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

void detect_collisions_from_candidates(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const Candidates& candidates,
    ConcurrentImpacts& impacts,
    TrajectoryType trajectory)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses + displacements);
        detect_collisions_from_candidates(
            V_t0, V_t1 - V_t0, bodies.m_edges, bodies.m_faces, candidates,
            impacts);
        return;
    }

    assert(trajectory == TrajectoryType::SCREWING);
    PROFILE_POINT("collisions_detection");
    NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

    PROFILE_START();
    PROFILE_START(NARROW_PHASE);

    auto detect_ev_collision = [&](const EdgeVertexCandidate& ev_candidate) {
        double toi, alpha;
        bool is_colliding = detect_edge_vertex_collisions_narrow_phase(
            bodies, poses, displacements, ev_candidate, toi, alpha, trajectory);
        if (is_colliding) {
            impacts.ev_impacts.emplace_back(
                toi, ev_candidate.edge_index, alpha, ev_candidate.vertex_index);
        }
    };

    auto detect_ee_collision = [&](const EdgeEdgeCandidate& ee_candidate) {
        double toi, edge0_alpha, edge1_alpha;
        bool is_colliding = detect_edge_edge_collisions_narrow_phase(
            bodies, poses, displacements, ee_candidate, toi, edge0_alpha,
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
            bodies, poses, displacements, fv_candidate, toi, u, v, trajectory);
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
    PROFILE_END();
}

// Determine if a single edge-vertext pair intersects.
bool detect_edge_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeVertexCandidate& candidate,
    double& toi,
    double& alpha,
    TrajectoryType trajectory)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        physics::Poses<double> poses_t1 = poses + displacements;
        Eigen::VectorX3d Vi_t0 =
            bodies.world_vertex(poses, bodies.m_edges(candidate.edge_index, 0));
        Eigen::VectorX3d Vi_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge_index, 0));

        Eigen::VectorX3d Vj_t0 =
            bodies.world_vertex(poses, bodies.m_edges(candidate.edge_index, 1));
        Eigen::VectorX3d Vj_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge_index, 1));

        Eigen::VectorX3d Vk_t0 =
            bodies.world_vertex(poses, candidate.vertex_index);
        Eigen::VectorX3d Vk_t1 =
            bodies.world_vertex(poses_t1, candidate.vertex_index);

        return detect_edge_vertex_collisions_narrow_phase(
            Vi_t0, Vj_t0, Vk_t0, Vi_t1 - Vi_t0, Vj_t1 - Vj_t0, Vk_t1 - Vk_t0,
            toi, alpha);
    }

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    bool is_colliding = compute_edge_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edge_id, toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA =
            poses[bodyA_id] + displacements[bodyA_id] * toi;
        physics::Pose<double> poseB =
            poses[bodyB_id] + displacements[bodyB_id] * toi;

        // Get the world vertex of the point at time t
        Eigen::VectorX3d v0 =
            bodies.m_rbs[bodyA_id].world_vertex(poseA, vertex_id);
        // Get the world vertex of the edge at time t
        Eigen::VectorX3d v1 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 0));
        Eigen::VectorX3d v2 = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edge_id, 1));

        // Compute the impact parameter along the edge
        alpha = geometry::project_point_to_line(v0, v1, (v2 - v1).eval());
        assert(
            alpha > -Constants::PARAMETER_ASSERTION_TOL
            && alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
    }
    return is_colliding;
}

bool detect_edge_edge_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const EdgeEdgeCandidate& candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    TrajectoryType trajectory)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        physics::Poses<double> poses_t1 = poses + displacements;

        Eigen::VectorX3d Vi_t0 = bodies.world_vertex(
            poses, bodies.m_edges(candidate.edge0_index, 0));
        Eigen::VectorX3d Vi_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge0_index, 0));

        Eigen::VectorX3d Vj_t0 = bodies.world_vertex(
            poses, bodies.m_edges(candidate.edge0_index, 1));
        Eigen::VectorX3d Vj_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge0_index, 1));

        Eigen::VectorX3d Vk_t0 = bodies.world_vertex(
            poses, bodies.m_edges(candidate.edge1_index, 0));
        Eigen::VectorX3d Vk_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge1_index, 0));

        Eigen::VectorX3d Vl_t0 = bodies.world_vertex(
            poses, bodies.m_edges(candidate.edge1_index, 1));
        Eigen::VectorX3d Vl_t1 = bodies.world_vertex(
            poses_t1, bodies.m_edges(candidate.edge1_index, 1));

        return detect_edge_edge_collisions_narrow_phase(
            Vi_t0, Vj_t0, Vk_t0, Vl_t0, Vi_t1 - Vi_t0, Vj_t1 - Vj_t0,
            Vk_t1 - Vk_t0, Vl_t1 - Vl_t0, toi, edge0_alpha, edge1_alpha);
    }

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);

    bool is_colliding = compute_edge_edge_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        edgeA_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], edgeB_id, toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA =
            poses[bodyA_id] + displacements[bodyA_id] * toi;
        physics::Pose<double> poseB =
            poses[bodyB_id] + displacements[bodyB_id] * toi;

        Eigen::Vector3d edgeA_vertex0_toi = bodies.m_rbs[bodyA_id].world_vertex(
            poseA, bodies.m_rbs[bodyA_id].edges(edgeA_id, 0));
        Eigen::Vector3d edgeA_vertex1_toi = bodies.m_rbs[bodyA_id].world_vertex(
            poseA, bodies.m_rbs[bodyA_id].edges(edgeA_id, 1));

        Eigen::Vector3d edgeB_vertex0_toi = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edgeB_id, 0));
        Eigen::Vector3d edgeB_vertex1_toi = bodies.m_rbs[bodyB_id].world_vertex(
            poseB, bodies.m_rbs[bodyB_id].edges(edgeB_id, 1));

        geometry::project_segment_to_segment(
            edgeA_vertex0_toi, edgeA_vertex1_toi, //
            edgeB_vertex0_toi, edgeB_vertex1_toi, //
            edge0_alpha, edge1_alpha);
        assert(
            edge0_alpha > -Constants::PARAMETER_ASSERTION_TOL
            && edge0_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
        assert(
            edge1_alpha > -Constants::PARAMETER_ASSERTION_TOL
            && edge1_alpha < 1 + Constants::PARAMETER_ASSERTION_TOL);
    }
    return is_colliding;
}

bool detect_face_vertex_collisions_narrow_phase(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses,
    const physics::Poses<double>& displacements,
    const FaceVertexCandidate& candidate,
    double& toi,
    double& u,
    double& v,
    TrajectoryType trajectory)
{
    if (trajectory == TrajectoryType::LINEARIZED) {
        physics::Poses<double> poses_t1 = poses + displacements;

        Eigen::VectorX3d Vi_t0 =
            bodies.world_vertex(poses, bodies.m_faces(candidate.face_index, 0));
        Eigen::VectorX3d Vi_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 0));

        Eigen::VectorX3d Vj_t0 =
            bodies.world_vertex(poses, bodies.m_faces(candidate.face_index, 1));
        Eigen::VectorX3d Vj_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 1));

        Eigen::VectorX3d Vk_t0 =
            bodies.world_vertex(poses, bodies.m_faces(candidate.face_index, 2));
        Eigen::VectorX3d Vk_t1 = bodies.world_vertex(
            poses_t1, bodies.m_faces(candidate.face_index, 2));

        Eigen::VectorX3d Vl_t0 =
            bodies.world_vertex(poses, candidate.vertex_index);
        Eigen::VectorX3d Vl_t1 =
            bodies.world_vertex(poses_t1, candidate.vertex_index);

        return detect_face_vertex_collisions_narrow_phase(
            Vi_t0, Vj_t0, Vk_t0, Vl_t0, Vi_t1 - Vi_t0, Vj_t1 - Vj_t0,
            Vk_t1 - Vk_t0, Vl_t1 - Vl_t0, toi, u, v);
    }

    assert(trajectory == TrajectoryType::SCREWING);
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);

    bool is_colliding = compute_face_vertex_time_of_impact(
        bodies.m_rbs[bodyA_id], poses[bodyA_id], displacements[bodyA_id],
        vertex_id, bodies.m_rbs[bodyB_id], poses[bodyB_id],
        displacements[bodyB_id], face_id, toi);
    if (is_colliding) {
        // Compute the poses at time toi
        physics::Pose<double> poseA =
            poses[bodyA_id] + displacements[bodyA_id] * toi;
        physics::Pose<double> poseB =
            poses[bodyB_id] + displacements[bodyB_id] * toi;

        // Get the world vertex of the point at time t
        Eigen::Vector3d vertex_toi =
            bodies.m_rbs[bodyA_id].world_vertex(poseA, vertex_id);
        // Get the world vertex of the face at time t
        Eigen::Vector3d triangle_vertex0_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 0));
        Eigen::Vector3d triangle_vertex1_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 1));
        Eigen::Vector3d triangle_vertex2_toi =
            bodies.m_rbs[bodyB_id].world_vertex(
                poseB, bodies.m_rbs[bodyB_id].faces(face_id, 2));

        // TODO: Consider moving this computation to an as needed basis
        double w;
        geometry::barycentric_coordinates(
            vertex_toi, triangle_vertex0_toi, triangle_vertex1_toi,
            triangle_vertex2_toi, u, v, w);
        assert(
            u + v + w < 1 + Constants::PARAMETER_ASSERTION_TOL
            && u + v + w > 1 - Constants::PARAMETER_ASSERTION_TOL);
    }
    return is_colliding;
}

} // namespace ccd
