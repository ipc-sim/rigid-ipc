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

namespace ipc::rigid {

typedef Pose<Interval> PoseI;

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

Eigen::Vector2d compute_edge_vertex_tolerance(
    const RigidBody& bodyA,       // Body of the vertex
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,             // In bodyA
    const RigidBody& bodyB,       // Body of the edge
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edge_id)               // In bodyB
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
        Constants::RIGID_CCD_TOI_TOL, // / dl,
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edge_id));
}

/// Find time-of-impact between two rigid bodies
bool compute_edge_vertex_time_of_impact(
    const RigidBody& bodyA,       // Body of the vertex
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,             // In bodyA
    const RigidBody& bodyB,       // Body of the edge
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edge_id,               // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 2);

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();

    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();

    const auto distance = [&](const VectorMax3I& params) {
        assert(params.size() == 2);
        return edge_vertex_aabb(
            bodyA, poseIA_t0, poseIA_t1, vertex_id, bodyB, poseIB_t0, poseIB_t1,
            edge_id, /*t=*/params(0), /*alpha=*/params(1));
    };

    Eigen::Vector2d tol = compute_edge_vertex_tolerance(
        bodyA, poseA_t0, poseA_t1, vertex_id, bodyB, poseB_t0, poseB_t1,
        edge_id);
    tol[0] = toi_tolerance;

    VectorMax3I x0 = Vector2I(Interval(0, earliest_toi), Interval(0, 1));
    VectorMax3I toi_interval;
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
    const RigidBody& bodyA,       // Body of the first edge
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,              // In bodyA
    const RigidBody& bodyB,       // Body of the second edge
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id)              // In bodyB
{

    std::cerr << fmt_eigen(bodyA.vertices.row(bodyA.edges(edgeA_id, 0)))
              << std::endl;
    std::cerr << fmt_eigen(bodyA.vertices.row(bodyA.edges(edgeA_id, 1)))
              << std::endl;
    std::cerr << fmt_eigen(poseA_t0.position) << std::endl;
    std::cerr << fmt_eigen(poseA_t0.rotation) << std::endl;
    std::cerr << fmt_eigen(poseA_t1.position) << std::endl;
    std::cerr << fmt_eigen(poseA_t1.rotation) << std::endl;
    std::cerr << fmt_eigen(bodyB.vertices.row(bodyB.edges(edgeB_id, 0)))
              << std::endl;
    std::cerr << fmt_eigen(bodyB.vertices.row(bodyB.edges(edgeB_id, 1)))
              << std::endl;
    std::cerr << fmt_eigen(poseB_t0.position) << std::endl;
    std::cerr << fmt_eigen(poseB_t0.rotation) << std::endl;
    std::cerr << fmt_eigen(poseB_t1.position) << std::endl;
    std::cerr << fmt_eigen(poseB_t1.rotation) << std::endl;
    std::cerr << std::endl;
}

Eigen::Vector3d compute_edge_edge_tolerance(
    const RigidBody& bodyA,       // Body of the first edge
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,              // In bodyA
    const RigidBody& bodyB,       // Body of the second edge
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id)              // In bodyB
{
    // TODO: Compute trajecotry length for interpolating rotation vectors.
    // double dl = ...

    return Eigen::Vector3d(
        Constants::RIGID_CCD_TOI_TOL, // / dl
        Constants::RIGID_CCD_LENGTH_TOL / bodyA.edge_length(edgeA_id),
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edgeB_id));
}

// Find time-of-impact between two rigid bodies
bool compute_edge_edge_time_of_impact(
    const RigidBody& bodyA,       // Body of the first edge
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t edgeA_id,              // In bodyA
    const RigidBody& bodyB,       // Body of the second edge
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t edgeB_id,              // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    assert(bodyA.dim() == 3 && bodyB.dim() == bodyA.dim());

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();
    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();
    const auto distance = [&](const VectorMax3I& params) {
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

    VectorMax3I toi_interval;
    VectorMax3I x0 =
        Vector3I(Interval(0, earliest_toi), Interval(0, 1), Interval(0, 1));
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

Eigen::Vector3d compute_face_vertex_tolerance(
    const RigidBody& bodyA,       // Body of the vertex
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,             // In bodyA
    const RigidBody& bodyB,       // Body of the triangle
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id)               // In bodyB
{
    // TODO: Compute trajecotry length for interpolating rotation vectors.
    // double dl = ...

    // (f1 - f0) * u + (f2 - f0) * v + f0
    // u interpolates edge (f0, f1) and v interpolates edge (f1, f2)
    size_t edge0_id = bodyB.mesh_selector.face_to_edge(face_id, 0);
    size_t edge1_id = bodyB.mesh_selector.face_to_edge(face_id, 1);

    return Eigen::Vector3d(
        // Constants::RIGID_CCD_TOI_TOL / dl,
        Constants::RIGID_CCD_TOI_TOL,
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edge0_id),
        Constants::RIGID_CCD_LENGTH_TOL / bodyB.edge_length(edge1_id));
}

// Find time-of-impact between two rigid bodies
bool compute_face_vertex_time_of_impact(
    const RigidBody& bodyA,       // Body of the vertex
    const Pose<double>& poseA_t0, // Pose of bodyA at t=0
    const Pose<double>& poseA_t1, // Pose of bodyA at t=1
    size_t vertex_id,             // In bodyA
    const RigidBody& bodyB,       // Body of the triangle
    const Pose<double>& poseB_t0, // Pose of bodyB at t=0
    const Pose<double>& poseB_t1, // Pose of bodyB at t=1
    size_t face_id,               // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    assert(bodyA.dim() == 3 && bodyA.dim() == bodyB.dim());

    const PoseI poseIA_t0 = poseA_t0.cast<Interval>();
    const PoseI poseIA_t1 = poseA_t1.cast<Interval>();
    const PoseI poseIB_t0 = poseB_t0.cast<Interval>();
    const PoseI poseIB_t1 = poseB_t1.cast<Interval>();

    const auto distance = [&](const VectorMax3I& params) {
        assert(params.size() == 3);
        return face_vertex_aabb(
            bodyA, poseIA_t0, poseIA_t1, vertex_id, //
            bodyB, poseIB_t0, poseIB_t1, face_id,   //
            /*t=*/params(0), /*u=*/params(1), /*v=*/params(2));
    };

    const auto is_domain_valid = [&](const VectorMax3I& params) {
        const Interval &t = params[0], &u = params[1], &v = params[2];
        // 0 ≤ t, u, v ≤ 1 is satisfied by the initial domain of the solve
        return overlap(u + v, Interval(0, 1));
    };

    Eigen::Vector3d tol = compute_face_vertex_tolerance(
        bodyA, poseA_t0, poseA_t1, vertex_id, //
        bodyB, poseB_t0, poseB_t1, face_id);
    tol[0] = toi_tolerance;

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    VectorMax3I toi_interval;
    VectorMax3I x0 =
        Vector3I(Interval(0, earliest_toi), Interval(0, 1), Interval(0, 1));
    bool is_impacting =
        interval_root_finder(distance, is_domain_valid, x0, tol, toi_interval);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    std::cout << "VF " << timer.getElapsedTime() << std::endl;
#endif

    // Return a conservative time-of-impact
    toi = is_impacting ? toi_interval(0).lower()
                       : std::numeric_limits<double>::infinity();

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ipc::rigid
