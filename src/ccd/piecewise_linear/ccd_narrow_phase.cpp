#include "ccd.hpp"

#include <igl/barycentric_coordinates.h>
#include <ipc/friction/closest_point.hpp>

#include <ccd/piecewise_linear/time_of_impact.hpp>
#include <constants.hpp>
#include <geometry/intersection.hpp>
#include <profiler.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

// Determine if a single edge-vertext pair intersects.
bool edge_vertex_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& candidate,
    double& toi,
    double& alpha,
    double earliest_toi)
{
    // TODO: Make this thread safe
    // PROFILE_POINT("detect_edge_vertex_collisions_narrow_phase");
    // PROFILE_START();

    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(candidate.edge_index, bodyB_id, edge_id);

    bool is_colliding = compute_piecewise_linear_edge_vertex_time_of_impact(
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

bool edge_edge_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& candidate,
    double& toi,
    double& edge0_alpha,
    double& edge1_alpha,
    double earliest_toi)
{
    // TODO: Make this thread safe
    // PROFILE_POINT("detect_edge_edge_collisions_narrow_phase");
    // PROFILE_START();

    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(candidate.edge1_index, bodyB_id, edgeB_id);

    bool is_colliding = compute_piecewise_linear_edge_edge_time_of_impact(
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

bool face_vertex_piecewise_linear_ccd(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& candidate,
    double& toi,
    double& u,
    double& v,
    double earliest_toi)
{
    // TODO: Make this thread safe
    // PROFILE_POINT("detect_face_vertex_collisions_narrow_phase");
    // PROFILE_START();

    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(candidate.face_index, bodyB_id, face_id);

    bool is_colliding = compute_piecewise_linear_face_vertex_time_of_impact(
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
