// Time-of-impact computation for rigid bodies with angular trajectories.
#include "rigid_body_time_of_impact.hpp"

#include <ccd/interval_root_finder.hpp>
#include <constants.hpp>
#include <geometry/distance.hpp>
#include <geometry/intersection.hpp>
#include <geometry/normal.hpp>
#include <logger.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

Eigen::Vector2d compute_edge_vertex_tolerance(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,       // Body of the edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& edge_id)                 // In bodyB
{
    double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();

    double omegaA_sqr = (poseA_t1.rotation - poseA_t0.rotation).squaredNorm();
    double omegaB_sqr = (poseB_t1.rotation - poseB_t0.rotation).squaredNorm();

    // Compute the maximum arc length of all the vertices
    double dl =
        sqrt(sA_sqr + omegaA_sqr * bodyA.vertices.row(vertex_id).squaredNorm());
    for (int i = 0; i < 2; i++) {
        double arc_len_sqr = omegaB_sqr
            * bodyB.vertices.row(bodyB.edges(edge_id, i)).squaredNorm();
        dl = std::max(dl, sqrt(sB_sqr + arc_len_sqr));
    }

    return Eigen::Vector2d(
        Constants::SCREWING_CCD_LENGTH_TOL / dl,
        Constants::SCREWING_CCD_LENGTH_TOL / bodyB.edge_length(edge_id));
}

/// Find time-of-impact between two rigid bodies
bool compute_edge_vertex_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,       // Body of the edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& edge_id,                 // In bodyB
    double& toi)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 2); // TODO: 3D

    typedef physics::Pose<Interval> PoseI;

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();

    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();

    const auto distance = [&](const Eigen::VectorXI& params) {
        assert(params.size() == 2);
        Interval t = params(0);
        Interval alpha = params(1);

        // Compute the poses at time t
        PoseI poseIA(
            (poseIA_t1.position - poseIA_t0.position) * t + poseIA_t0.position,
            (poseIA_t1.rotation - poseIA_t0.rotation) * t + poseIA_t0.rotation);
        PoseI poseIB(
            (poseIB_t1.position - poseIB_t0.position) * t + poseIB_t0.position,
            (poseIB_t1.rotation - poseIB_t0.rotation) * t + poseIB_t0.rotation);

        // Get the world vertex of the edges at time t
        Eigen::Vector2I vertex = bodyA.world_vertex(poseIA, vertex_id);

        // Get the world vertex of the edge at time t
        Eigen::Vector2I edge_vertex0 =
            bodyB.world_vertex(poseIB, bodyB.edges(edge_id, 0));
        Eigen::Vector2I edge_vertex1 =
            bodyB.world_vertex(poseIB, bodyB.edges(edge_id, 1));
        Eigen::Vector2I edge_vertex =
            (edge_vertex1 - edge_vertex0) * alpha + edge_vertex0;

        return (vertex - edge_vertex).eval();
    };

    Eigen::Vector2d tol = compute_edge_vertex_tolerance(
        bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
        edge_id);

    Eigen::VectorXI toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::VectorXI::Constant(2, Interval(0, 1)), tol,
        toi_interval);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

Eigen::Vector3d compute_edge_edge_tolerance(
    const physics::RigidBody& bodyA,       // Body of the first edge
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    double omegaA,
    const Eigen::Matrix3d& RA,
    const size_t& edgeA_id,                // In bodyA
    const physics::RigidBody& bodyB,       // Body of the second edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    double omegaB,
    const Eigen::Matrix3d& RB,
    const size_t& edgeB_id) // In bodyB
{
    double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();

    double dl = -std::numeric_limits<double>::infinity();
    Eigen::Vector3d v;
    double radius_sqr;
    for (int i = 0; i < 2; i++) {
        v = bodyA.vertices.row(bodyA.edges(edgeA_id, i));
        radius_sqr = (RA * v).head<2>().squaredNorm();
        dl = std::max(dl, sqrt(sA_sqr + omegaA * omegaA * radius_sqr));

        v = bodyB.vertices.row(bodyB.edges(edgeB_id, i));
        radius_sqr = (RB * v).head<2>().squaredNorm();
        dl = std::max(dl, sqrt(sB_sqr + omegaB * omegaB * radius_sqr));
    }

    return Eigen::Vector3d(
        Constants::SCREWING_CCD_LENGTH_TOL / dl,
        Constants::SCREWING_CCD_LENGTH_TOL / bodyA.edge_length(edgeA_id),
        Constants::SCREWING_CCD_LENGTH_TOL / bodyB.edge_length(edgeB_id));
}

// Find time-of-impact between two rigid bodies
bool compute_edge_edge_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the first edge
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& edgeA_id,                // In bodyA
    const physics::RigidBody& bodyB,       // Body of the second edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& edgeB_id,                // In bodyB
    double& toi)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

    const Eigen::Vector3I posA_t0 = poseA_t0.position.cast<Interval>();
    const Eigen::Vector3I posA_t1 = poseA_t1.position.cast<Interval>();
    const Eigen::Vector3I posB_t0 = poseB_t0.position.cast<Interval>();
    const Eigen::Vector3I posB_t1 = poseB_t1.position.cast<Interval>();

    Eigen::Matrix3d RA_t0, PA;
    double omegaA;
    decompose_to_z_screwing(poseA_t0, poseA_t1, RA_t0, PA, omegaA);
    Eigen::Matrix3I RA_t0_I = RA_t0.cast<Interval>();
    Eigen::Matrix3I PA_I = PA.cast<Interval>();

    Eigen::Matrix3d RB_t0, PB;
    double omegaB;
    decompose_to_z_screwing(poseB_t0, poseB_t1, RB_t0, PB, omegaB);
    Eigen::Matrix3I RB_t0_I = RB_t0.cast<Interval>();
    Eigen::Matrix3I PB_I = PB.cast<Interval>();

    const auto distance = [&](const Eigen::VectorXI& params) {
        assert(params.size() == 3);
        Interval t = params(0);
        Interval edgeA_alpha = params(1);
        Interval edgeB_alpha = params(2);

        // Compute the poses at time t
        Eigen::VectorX3<Interval> posA_t = (posA_t1 - posA_t0) * t + posA_t0;
        Eigen::MatrixXX3<Interval> RA_t =
            PA_I.transpose() * rotate_around_z(omegaA * t) * PA_I * RA_t0_I;

        Eigen::VectorX3<Interval> posB_t = (posB_t1 - posB_t0) * t + posB_t0;
        Eigen::MatrixXX3<Interval> RB_t =
            PB_I.transpose() * rotate_around_z(omegaB * t) * PB_I * RB_t0_I;

        // Get the world vertex of the edges at time t
        Eigen::Vector3I edgeA_vertex0 =
            bodyA.world_vertex(RA_t, posA_t, bodyA.edges(edgeA_id, 0));
        Eigen::Vector3I edgeA_vertex1 =
            bodyA.world_vertex(RA_t, posA_t, bodyA.edges(edgeA_id, 1));
        Eigen::Vector3I edgeA_vertex =
            (edgeA_vertex1 - edgeA_vertex0) * edgeA_alpha + edgeA_vertex0;

        Eigen::Vector3I edgeB_vertex0 =
            bodyB.world_vertex(RB_t, posB_t, bodyB.edges(edgeB_id, 0));
        Eigen::Vector3I edgeB_vertex1 =
            bodyB.world_vertex(RB_t, posB_t, bodyB.edges(edgeB_id, 1));
        Eigen::Vector3I edgeB_vertex =
            (edgeB_vertex1 - edgeB_vertex0) * edgeB_alpha + edgeB_vertex0;

        return (edgeB_vertex - edgeA_vertex).eval();
    };

    Eigen::Vector3d tol = compute_edge_edge_tolerance(
        bodyA, poseA_t0, poseA_t1, omegaA, PA * RA_t0, edgeA_id, //
        bodyB, poseB_t0, poseB_t1, omegaB, PB * RB_t0, edgeB_id);

    Eigen::VectorXI toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::VectorXI::Constant(3, Interval(0, 1)), tol,
        toi_interval);
    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Face-Vertex

double compute_face_vertex_tolerance(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    double omegaA,
    const Eigen::Matrix3d& RA,
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,       // Body of the triangle
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    double omegaB,
    const Eigen::Matrix3d& RB,
    const size_t& face_id) // In bodyB
{
    double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();

    Eigen::Vector3d v = bodyA.vertices.row(vertex_id);
    double radius_sqr = (RA * v).head<2>().squaredNorm();
    double dl = sqrt(sA_sqr + omegaA * omegaA * v.squaredNorm());

    for (int i = 0; i < 3; i++) {
        v = bodyB.vertices.row(bodyB.faces(face_id, i));
        radius_sqr = (RB * v).head<2>().squaredNorm();
        dl = std::max(dl, sqrt(sB_sqr + omegaB * omegaB * radius_sqr));
    }

    return Constants::SCREWING_CCD_LENGTH_TOL / dl;
}

// Find time-of-impact between two rigid bodies
bool compute_face_vertex_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    const size_t& vertex_id,               // In bodyA
    const physics::RigidBody& bodyB,       // Body of the triangle
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    const size_t& face_id,                 // In bodyB
    double& toi)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

    const Eigen::Vector3I posA_t0 = poseA_t0.position.cast<Interval>();
    const Eigen::Vector3I posA_t1 = poseA_t1.position.cast<Interval>();
    const Eigen::Vector3I posB_t0 = poseB_t0.position.cast<Interval>();
    const Eigen::Vector3I posB_t1 = poseB_t1.position.cast<Interval>();

    Eigen::Matrix3d RA_t0, PA;
    double omegaA;
    decompose_to_z_screwing(poseA_t0, poseA_t1, RA_t0, PA, omegaA);
    Eigen::Matrix3I RA_t0_I = RA_t0.cast<Interval>();
    Eigen::Matrix3I PA_I = PA.cast<Interval>();

    Eigen::Matrix3d RB_t0, PB;
    double omegaB;
    decompose_to_z_screwing(poseB_t0, poseB_t1, RB_t0, PB, omegaB);
    Eigen::Matrix3I RB_t0_I = RB_t0.cast<Interval>();
    Eigen::Matrix3I PB_I = PB.cast<Interval>();

    const auto vertex_positions = [&](const Interval& t,
                                      Eigen::Vector3I& vertex,
                                      Eigen::Vector3I& face_vertex0,
                                      Eigen::Vector3I& face_vertex1,
                                      Eigen::Vector3I& face_vertex2) {
        // Compute the poses at time t
        Eigen::VectorX3<Interval> posA_t = (posA_t1 - posA_t0) * t + posA_t0;
        Eigen::MatrixXX3<Interval> RA_t =
            PA_I.transpose() * rotate_around_z(omegaA * t) * PA_I * RA_t0_I;

        Eigen::VectorX3<Interval> posB_t = (posB_t1 - posB_t0) * t + posB_t0;
        Eigen::MatrixXX3<Interval> RB_t =
            PB_I.transpose() * rotate_around_z(omegaB * t) * PB_I * RB_t0_I;

        // Get the world vertex of the point at time t
        vertex = bodyA.world_vertex(RA_t, posA_t, vertex_id);
        // Get the world vertex of the edge at time t
        face_vertex0 =
            bodyB.world_vertex(RB_t, posB_t, bodyB.faces(face_id, 0));
        face_vertex1 =
            bodyB.world_vertex(RB_t, posB_t, bodyB.faces(face_id, 1));
        face_vertex2 =
            bodyB.world_vertex(RB_t, posB_t, bodyB.faces(face_id, 2));
    };

    const auto distance = [&](const Interval& t) {
        // Get the world vertex and face of the point at time t
        Eigen::Vector3I vertex;
        Eigen::Vector3I face_vertex0, face_vertex1, face_vertex2;
        vertex_positions(t, vertex, face_vertex0, face_vertex1, face_vertex2);
        return geometry::point_plane_signed_distance(
            vertex, face_vertex0, face_vertex1, face_vertex2);
    };

    const auto is_point_inside_triangle = [&](const Interval& t) {
        // Get the world vertex and face of the point at time t
        Eigen::Vector3I vertex;
        Eigen::Vector3I face_vertex0, face_vertex1, face_vertex2;
        vertex_positions(t, vertex, face_vertex0, face_vertex1, face_vertex2);
        return geometry::is_point_inside_triangle(
            vertex, face_vertex0, face_vertex1, face_vertex2);
    };

    double tol = compute_face_vertex_tolerance(
        bodyA, poseA_t0, poseA_t1, omegaA, PA * RA_t0, vertex_id, //
        bodyB, poseB_t0, poseB_t1, omegaB, PB * RB_t0, face_id);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_point_inside_triangle, Interval(0, 1), tol, toi_interval);
    // Return a conservative time-of-impact
    toi = toi_interval.lower();
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ccd
