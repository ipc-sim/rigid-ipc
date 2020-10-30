// Time-of-impact computation for rigid bodies with angular trajectories.
#include "time_of_impact.hpp"

// #define TIME_CCD_QUERIES
#ifdef TIME_CCD_QUERIES
#include <igl/Timer.h>
#endif

#include <ccd/rigid/rigid_trajectory_aabb.hpp>
#include <geometry/distance.hpp>
#include <geometry/intersection.hpp>
#include <interval/interval_root_finder.hpp>
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
    size_t vertex_id,                      // In bodyA
    const physics::RigidBody& bodyB,       // Body of the edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edge_id)                        // In bodyB
{
    // double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    // double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();
    //
    // double omegaA_sqr = (poseA_t1.rotation -
    // poseA_t0.rotation).squaredNorm(); double omegaB_sqr = (poseB_t1.rotation
    // - poseB_t0.rotation).squaredNorm();
    //
    // // Compute the maximum arc length of all the vertices
    // double dl =
    //     sqrt(sA_sqr + omegaA_sqr *
    //     bodyA.vertices.row(vertex_id).squaredNorm());
    // for (int i = 0; i < 2; i++) {
    //     double arc_len_sqr = omegaB_sqr
    //         * bodyB.vertices.row(bodyB.edges(edge_id, i)).squaredNorm();
    //     dl = std::max(dl, sqrt(sB_sqr + arc_len_sqr));
    // }

    return Eigen::Vector2d(
        // 1e-6 / dl,
        // 1e-6 / bodyB.edge_length(edge_id));
        Constants::RIGID_CCD_TOI_TOL,
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edge_id));
}

/// Find time-of-impact between two rigid bodies
bool compute_edge_vertex_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,                      // In bodyA
    const physics::RigidBody& bodyB,       // Body of the edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edge_id,                        // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 2); // TODO: 3D

    typedef physics::Pose<Interval> PoseI;

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();

    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();

    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 2);
        return edge_vertex_aabb(
            bodyA, poseIA_t0, poseIA_t1, vertex_id, bodyB, poseIB_t0, poseIB_t1,
            edge_id, /*t=*/params(0), /*alpha=*/params(1));
    };

    Eigen::Vector2d tol = compute_edge_vertex_tolerance(
        bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
        edge_id);
    tol[0] = toi_tolerance;

    Eigen::VectorX3I x0 =
        Eigen::Vector2I(Interval(0, earliest_toi), Interval(0, 1));
    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(distance, x0, tol, toi_interval);

    // Return a conservative time-of-impact
    toi = is_impacting ? toi_interval(0).lower()
                       : std::numeric_limits<double>::infinity();

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

void print_EE_query(
    const physics::RigidBody& bodyA,       // Body of the first edge
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,                       // In bodyA
    const physics::RigidBody& bodyB,       // Body of the second edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id)                       // In bodyB
{

    std::cerr << logger::fmt_eigen(bodyA.vertices.row(bodyA.edges(edgeA_id, 0)))
              << std::endl;
    std::cerr << logger::fmt_eigen(bodyA.vertices.row(bodyA.edges(edgeA_id, 1)))
              << std::endl;
    std::cerr << logger::fmt_eigen(poseA_t0.position) << std::endl;
    std::cerr << logger::fmt_eigen(poseA_t0.rotation) << std::endl;
    std::cerr << logger::fmt_eigen(poseA_t1.position) << std::endl;
    std::cerr << logger::fmt_eigen(poseA_t1.rotation) << std::endl;
    std::cerr << logger::fmt_eigen(bodyB.vertices.row(bodyB.edges(edgeB_id, 0)))
              << std::endl;
    std::cerr << logger::fmt_eigen(bodyB.vertices.row(bodyB.edges(edgeB_id, 1)))
              << std::endl;
    std::cerr << logger::fmt_eigen(poseB_t0.position) << std::endl;
    std::cerr << logger::fmt_eigen(poseB_t0.rotation) << std::endl;
    std::cerr << logger::fmt_eigen(poseB_t1.position) << std::endl;
    std::cerr << logger::fmt_eigen(poseB_t1.rotation) << std::endl;
    std::cerr << std::endl;
}

Eigen::Vector3d compute_edge_edge_tolerance(
    const physics::RigidBody& bodyA,       // Body of the first edge
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,                       // In bodyA
    const physics::RigidBody& bodyB,       // Body of the second edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id)                       // In bodyB
{
    // TODO: This is not exactly the arc length of the trajectory
    // Eigen::Matrix3d RA_t0, PA;
    // double omegaA;
    // decompose_to_z_screwing(poseA_t0, poseA_t1, RA_t0, PA, omegaA);
    // Eigen::Matrix3d RA = PA * RA_t0;
    //
    // Eigen::Matrix3d RB_t0, PB;
    // double omegaB;
    // decompose_to_z_screwing(poseB_t0, poseB_t1, RB_t0, PB, omegaB);
    // Eigen::Matrix3d RB = PB * RB_t0;
    //
    // double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    // double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();
    //
    // double dl = -std::numeric_limits<double>::infinity();
    // Eigen::Vector3d v;
    // double radius_sqr;
    // for (int i = 0; i < 2; i++) {
    //     v = bodyA.vertices.row(bodyA.edges(edgeA_id, i));
    //     radius_sqr = (RA * v).head<2>().squaredNorm();
    //     dl = std::max(dl, sqrt(sA_sqr + omegaA * omegaA * radius_sqr));
    //
    //     v = bodyB.vertices.row(bodyB.edges(edgeB_id, i));
    //     radius_sqr = (RB * v).head<2>().squaredNorm();
    //     dl = std::max(dl, sqrt(sB_sqr + omegaB * omegaB * radius_sqr));
    // }

    return Eigen::Vector3d(
        // 1e-6 / dl,
        // 1e-6 / bodyA.edge_length(edgeA_id),
        // 1e-6 / bodyB.edge_length(edgeB_id));
        Constants::RIGID_CCD_TOI_TOL,
        Constants::RIGID_CCD_LENGTH_TOL / bodyA.edge_length(edgeA_id),
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edgeB_id));
}

// Find time-of-impact between two rigid bodies
bool compute_edge_edge_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the first edge
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,                       // In bodyA
    const physics::RigidBody& bodyB,       // Body of the second edge
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id,                       // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    assert(bodyA.dim() == 3 && bodyB.dim() == bodyA.dim());

    const physics::Pose<Interval> poseIA_t0 = poseA_t0.cast<Interval>();
    const physics::Pose<Interval> poseIA_t1 = poseA_t1.cast<Interval>();
    const physics::Pose<Interval> poseIB_t0 = poseB_t0.cast<Interval>();
    const physics::Pose<Interval> poseIB_t1 = poseB_t1.cast<Interval>();
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 3);
        return edge_edge_aabb(
            bodyA, poseIA_t0, poseIA_t1, edgeA_id, bodyB, poseIB_t0, poseIB_t1,
            edgeB_id, /*t=*/params(0), /*alpha=*/params(1), /*beta=*/params(2));
    };

    Eigen::Vector3d tol = compute_edge_edge_tolerance(
        bodyA, poseA_t0, poseA_t1, edgeA_id, //
        bodyB, poseB_t0, poseB_t1, edgeB_id);
    tol[0] = toi_tolerance;

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    Eigen::VectorX3I toi_interval;
    Eigen::VectorX3I x0 = Eigen::Vector3I(
        Interval(0, earliest_toi), Interval(0, 1), Interval(0, 1));
    bool is_impacting = interval_root_finder(distance, x0, tol, toi_interval);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    std::cout << "EE " << timer.getElapsedTime() << std::endl;
    if (timer.getElapsedTime() >= 60) {
        std::cerr << "EE (" << timer.getElapsedTime() << "s)" << std::endl;
        print_EE_query(
            bodyA, poseA_t0, poseA_t1, edgeA_id, //
            bodyB, poseB_t0, poseB_t1, edgeB_id);
    }
#endif

    // Return a conservative time-of-impact
    toi = is_impacting ? toi_interval(0).lower()
                       : std::numeric_limits<double>::infinity();

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
    size_t vertex_id,                      // In bodyA
    const physics::RigidBody& bodyB,       // Body of the triangle
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id)                        // In bodyB
{
    // TODO: This is not exactly the arc length of the trajectory
    // Eigen::Matrix3d RA_t0, PA;
    // double omegaA;
    // decompose_to_z_screwing(poseA_t0, poseA_t1, RA_t0, PA, omegaA);
    // Eigen::Matrix3d RA = PA * RA_t0;
    //
    // Eigen::Matrix3d RB_t0, PB;
    // double omegaB;
    // decompose_to_z_screwing(poseB_t0, poseB_t1, RB_t0, PB, omegaB);
    // Eigen::Matrix3d RB = PB * RB_t0;
    //
    // double sA_sqr = (poseA_t1.position - poseA_t0.position).squaredNorm();
    // double sB_sqr = (poseB_t1.position - poseB_t0.position).squaredNorm();
    //
    // Eigen::Vector3d v = bodyA.vertices.row(vertex_id);
    // double radius_sqr = (RA * v).head<2>().squaredNorm();
    // double dl = sqrt(sA_sqr + omegaA * omegaA * v.squaredNorm());
    //
    // for (int i = 0; i < 3; i++) {
    //     v = bodyB.vertices.row(bodyB.faces(face_id, i));
    //     radius_sqr = (RB * v).head<2>().squaredNorm();
    //     dl = std::max(dl, sqrt(sB_sqr + omegaB * omegaB * radius_sqr));
    // }

    // return 1e-6 / dl;
    return Constants::RIGID_CCD_TOI_TOL;
}

// Find time-of-impact between two rigid bodies
bool compute_face_vertex_time_of_impact(
    const physics::RigidBody& bodyA,       // Body of the vertex
    const physics::Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const physics::Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,                      // In bodyA
    const physics::RigidBody& bodyB,       // Body of the triangle
    const physics::Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const physics::Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id,                        // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

    typedef physics::Pose<Interval> PoseI;

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();

    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();

    const auto vertex_positions =
        [&](const Interval& t, Eigen::Vector3I& vertex,
            Eigen::Vector3I& face_vertex0, Eigen::Vector3I& face_vertex1,
            Eigen::Vector3I& face_vertex2) {
            // Compute the poses at time t
            PoseI poseIA = PoseI::interpolate(poseIA_t0, poseIA_t1, t);
            PoseI poseIB = PoseI::interpolate(poseIB_t0, poseIB_t1, t);

            // Get the world vertex of the point at time t
            vertex = bodyA.world_vertex(poseIA, vertex_id);
            // Get the world vertex of the edge at time t
            face_vertex0 = bodyB.world_vertex(poseIB, bodyB.faces(face_id, 0));
            face_vertex1 = bodyB.world_vertex(poseIB, bodyB.faces(face_id, 1));
            face_vertex2 = bodyB.world_vertex(poseIB, bodyB.faces(face_id, 2));
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

    // double tol = compute_face_vertex_tolerance(
    //     bodyA, poseA_t0, poseA_t1, vertex_id, //
    //     bodyB, poseB_t0, poseB_t1, face_id);
    double tol = toi_tolerance;

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_point_inside_triangle, Interval(0, earliest_toi), tol,
        toi_interval);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    std::cout << "VF " << timer.getElapsedTime() << std::endl;
#endif

    // Return a conservative time-of-impact
    toi = is_impacting ? toi_interval.lower()
                       : std::numeric_limits<double>::infinity();

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ccd
