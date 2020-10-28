#include "ccd.hpp"

#include <ccd/linear/ccd.hpp>
#include <ccd/piecewise_linear/ccd.hpp>
#include <ccd/rigid/ccd.hpp>

namespace ccd {

void detect_collisions(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ConcurrentImpacts& impacts,
    DetectionMethod method,
    TrajectoryType trajectory)
{
    assert(bodies.num_bodies() == poses_t0.size());
    assert(poses_t0.size() == poses_t1.size());

    // Do the broad phase by detecting candidate impacts
    ipc::Candidates candidates;
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
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const int collision_types,
    ipc::Candidates& candidates,
    DetectionMethod method,
    TrajectoryType trajectory,
    const double inflation_radius)
{
    if (bodies.m_rbs.size() <= 1) {
        return;
    }

    // TODO: Add a quick capsule test for intersecting rigid bodies
    // def capsule(body, pose, displacement):
    //     return {"radius": body.r_max, "start": pose.position, "end":
    //     pose.position + displacement.position}
    //
    // def distance(capsule1, capsule2):
    //     return segment_segment_distance(capsule1.start, capsule1.end,
    //     capsule2.start, capsule2.end) - capsule1.radius - capsule2.radius
    //
    // capsules = [capsule(body, pose, displacement) for body in bodies]
    //
    // intersecting_capsules = {}
    // # O(n^2) segment_segment_distance calls
    // for i in range(len(capsules)):
    //     for j in range(i, len(capsules)):
    //         if distance(capsules[i], capsules[j]) <= tol:
    //             intersecting_capsules.add(i)
    //             intersecting_capsules.add(j)
    //
    // detect_collisions(bodies[intersecting_capsules], ...)

    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses_t0);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses_t1);
        detect_collision_candidates(
            V_t0, V_t1, bodies.m_edges, bodies.m_faces, bodies.group_ids(),
            collision_types, candidates, method, inflation_radius);
        break;
    }
    case TrajectoryType::PIECEWISE_LINEAR:
    case TrajectoryType::RIGID:
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
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::Candidates& candidates,
    ConcurrentImpacts& impacts,
    TrajectoryType trajectory)
{
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
        Eigen::MatrixXd V_t0 = bodies.world_vertices(poses_t0);
        Eigen::MatrixXd V_t1 = bodies.world_vertices(poses_t1);
        detect_collisions_from_candidates(
            V_t0, V_t1, bodies.m_edges, bodies.m_faces, candidates, impacts);
        break;
    }
    case TrajectoryType::PIECEWISE_LINEAR:
        throw NotImplementedError(
            "detect_collisions_from_candidates() is not implemented for "
            "piecewise linear trajectories!");
        break;
    case TrajectoryType::RIGID:
        detect_rigid_collisions_from_candidates(
            bodies, poses_t0, poses_t1, candidates, impacts);
        break;
    }
}

// Determine if a single edge-vertext pair intersects.
bool edge_vertex_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& candidate,
    double& toi,
    double& alpha,
    TrajectoryType trajectory,
    double earliest_toi)
{
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
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

    case TrajectoryType::PIECEWISE_LINEAR:
        return edge_vertex_piecewise_linear_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, alpha, earliest_toi);

    case TrajectoryType::RIGID:
        return edge_vertex_rigid_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, alpha, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

bool edge_edge_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    TrajectoryType trajectory,
    double earliest_toi)
{
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
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

    case TrajectoryType::PIECEWISE_LINEAR:
        return edge_edge_piecewise_linear_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, edge0_alpha,
            edge1_alpha, earliest_toi);

    case TrajectoryType::RIGID:
        return edge_edge_rigid_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, edge0_alpha,
            edge1_alpha, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

bool face_vertex_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& candidate,
    double& toi,
    double& u,
    double& v,
    TrajectoryType trajectory,
    double earliest_toi)
{
    switch (trajectory) {
    case TrajectoryType::LINEAR: {
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

    case TrajectoryType::PIECEWISE_LINEAR:
        return face_vertex_piecewise_linear_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, u, v, earliest_toi);

    case TrajectoryType::RIGID:
        return face_vertex_rigid_ccd(
            bodies, poses_t0, poses_t1, candidate, toi, u, v, earliest_toi);

    default:
        throw "Invalid trajectory type";
    }
}

} // namespace ccd
