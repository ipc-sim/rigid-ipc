// Time-of-impact computation for rigid bodies with angular trajectories.
#include "time_of_impact.hpp"

#include <stack>

#include <tight_inclusion/inclusion_ccd.hpp>

// #define TIME_CCD_QUERIES
#ifdef TIME_CCD_QUERIES
#include <igl/Timer.h>
#endif

// #define USE_FIXED_PIECES

#include <ipc/ipc.hpp>

#include <ccd/linear/edge_vertex_ccd.hpp>
#include <ccd/rigid/rigid_trajectory_aabb.hpp>
#include <interval/interval.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

typedef physics::Pose<Interval> PoseI;
typedef physics::Pose<double> Pose;

static const int NUM_PIECES = 100;

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

/// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_vertex_time_of_impact(
    const physics::RigidBody& bodyA, // Body of the vertex
    const Pose& poseA_t0,            // Pose of bodyA at t=0
    const Pose& poseA_t1,            // Pose of bodyA at t=1
    size_t vertex_id,                // In bodyA
    const physics::RigidBody& bodyB, // Body of the edge
    const Pose& poseB_t0,            // Pose of bodyB at t=0
    const Pose& poseB_t1,            // Pose of bodyB at t=1
    size_t edge_id,                  // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 2);

    spdlog::warn("piecewise linear edge vertex CCD is not conservative");

    bool is_impacting = false;
    Pose poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    for (int i = 1; i <= NUM_PIECES; i++) {
        double ti0 = (i - 1) / double(NUM_PIECES) * earliest_toi;
        double ti1 = i / double(NUM_PIECES) * earliest_toi;
        Pose poseA_ti1 = Pose::interpolate(poseA_t0, poseA_t1, ti1);
        Pose poseB_ti1 = Pose::interpolate(poseB_t0, poseB_t1, ti1);

        Eigen::Vector2d v_t0 = bodyA.world_vertex(poseA_ti0, vertex_id);
        Eigen::Vector2d v_t1 = bodyA.world_vertex(poseA_ti1, vertex_id);
        Eigen::Vector2d e0_t0 =
            bodyA.world_vertex(poseB_ti0, bodyB.edges(edge_id, 0));
        Eigen::Vector2d e0_t1 =
            bodyA.world_vertex(poseB_ti1, bodyB.edges(edge_id, 0));
        Eigen::Vector2d e1_t0 =
            bodyA.world_vertex(poseB_ti0, bodyB.edges(edge_id, 1));
        Eigen::Vector2d e1_t1 =
            bodyA.world_vertex(poseB_ti1, bodyB.edges(edge_id, 1));

        is_impacting = autodiff::compute_edge_vertex_time_of_impact(
            e0_t0, e1_t1, v_t0, (e0_t1 - e0_t0).eval(), (e1_t1 - e1_t0).eval(),
            (v_t1 - v_t0).eval(), toi);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            break;
        }

        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_edge_time_of_impact(
    const physics::RigidBody& bodyA, // Body of the first edge
    const Pose& poseA_t0,            // Pose of bodyA at t=0
    const Pose& poseA_t1,            // Pose of bodyA at t=1
    size_t edgeA_id,                 // In bodyA
    const physics::RigidBody& bodyB, // Body of the second edge
    const Pose& poseB_t0,            // Pose of bodyB at t=0
    const Pose& poseB_t1,            // Pose of bodyB at t=1
    size_t edgeB_id,                 // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

    long ea0i = bodyA.edges(edgeA_id, 0);
    long ea1i = bodyA.edges(edgeA_id, 1);
    long eb0i = bodyB.edges(edgeB_id, 0);
    long eb1i = bodyB.edges(edgeB_id, 1);

    double distance_t0 = sqrt(ipc::edge_edge_distance(
        bodyA.world_vertex(poseA_t0, ea0i), bodyA.world_vertex(poseA_t0, ea1i),
        bodyB.world_vertex(poseB_t0, eb0i),
        bodyB.world_vertex(poseB_t0, eb1i)));

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    Pose poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = NUM_PIECES; i > 0; i--) {
        ts.push(i / double(NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        Pose poseA_ti1 = Pose::interpolate(poseA_t0, poseA_t1, ti1);
        Pose poseB_ti1 = Pose::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(ipc::edge_edge_distance(
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 1))));

        if (distance_ti0 < 0.2 * distance_t0 && ti0 >= 1e-6) {
            toi = ti0;
            is_impacting = true;
            break;
        }

        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval t(ti0, ti1);

        PoseI poseIA_t0 = poseA_t0.cast<Interval>();
        PoseI poseIA_t1 = poseA_t1.cast<Interval>();
        PoseI poseIB_t0 = poseB_t0.cast<Interval>();
        PoseI poseIB_t1 = poseB_t1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_t0, poseIA_t1, t);
        PoseI poseIB = PoseI::interpolate(poseIB_t0, poseIB_t1, t);

        Eigen::MatrixXX3I RA = poseIA.construct_rotation_matrix();
        Eigen::MatrixXX3I RB = poseIB.construct_rotation_matrix();

        using Vector3I = Eigen::Vector3I;

        Vector3I ea0_t0 = bodyA.world_vertex(poseA_t0, ea0i).cast<Interval>();
        Vector3I ea0_t1 = bodyA.world_vertex(poseA_t1, ea0i).cast<Interval>();
        Vector3I ea0 = bodyA.world_vertex(RA, poseIA.position, ea0i);
        Interval d = (ea0 - ((ea0_t1 - ea0_t0) * t + ea0_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I ea1_t0 = bodyA.world_vertex(poseA_t0, ea1i).cast<Interval>();
        Vector3I ea1_t1 = bodyA.world_vertex(poseA_t1, ea1i).cast<Interval>();
        Vector3I ea1 = bodyA.world_vertex(RA, poseIA.position, ea1i);
        d = (ea1 - ((ea1_t1 - ea1_t0) * t + ea1_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I eb0_t0 = bodyB.world_vertex(poseB_t0, eb0i).cast<Interval>();
        Vector3I eb0_t1 = bodyB.world_vertex(poseB_t1, eb0i).cast<Interval>();
        Vector3I eb0 = bodyB.world_vertex(RB, poseIB.position, eb0i);
        d = (eb0 - ((eb0_t1 - eb0_t0) * t + eb0_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I eb1_t0 = bodyB.world_vertex(poseB_t0, eb1i).cast<Interval>();
        Vector3I eb1_t1 = bodyB.world_vertex(poseB_t1, eb1i).cast<Interval>();
        Vector3I eb1 = bodyB.world_vertex(RB, poseIB.position, eb1i);
        d = (eb1 - ((eb1_t1 - eb1_t0) * t + eb1_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        if (min_distance >= 0.5 * distance_ti0) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        // spdlog::critical("min_distance={:g}", min_distance);

        double output_tolerance;
        // 0: normal ccd method which only checks t = [0,1]
        // 1: ccd with max_itr and t=[0, t_max]
        const int CCD_TYPE = 1;
        is_impacting = inclusion_ccd::edgeEdgeCCD_double(
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 1)),
            bodyA.world_vertex(poseA_ti1, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti1, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti1, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti1, bodyB.edges(edgeB_id, 1)),
            { { -1, -1, -1 } }, // rounding error
            min_distance,       // minimum separation distance
            toi,                // time of impact
            1e-6,               // delta
            1.0,                // Maximum time to check
            1e6,                // Maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            break;
        }

        ts.pop();
        ti0 = ti1;
        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }
    // spdlog::trace("ee_ccd_num_subdivision={:d}", num_subdivisions);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    fmt::print("EE {:.16g}\n", timer.getElapsedTime());
#endif

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Face-Vertex

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_face_vertex_time_of_impact(
    const physics::RigidBody& bodyA, // Body of the vertex
    const Pose& poseA_t0,            // Pose of bodyA at t=0
    const Pose& poseA_t1,            // Pose of bodyA at t=1
    size_t vertex_id,                // In bodyA
    const physics::RigidBody& bodyB, // Body of the triangle
    const Pose& poseB_t0,            // Pose of bodyB at t=0
    const Pose& poseB_t1,            // Pose of bodyB at t=1
    size_t face_id,                  // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

    typedef physics::Pose<Interval> PoseI;

    const auto& vi = vertex_id;
    const auto& f0i = bodyB.faces(face_id, 0);
    const auto& f1i = bodyB.faces(face_id, 1);
    const auto& f2i = bodyB.faces(face_id, 2);

    double distance_t0 = sqrt(ipc::point_triangle_distance(
        bodyA.world_vertex(poseA_t0, vi), bodyB.world_vertex(poseB_t0, f0i),
        bodyB.world_vertex(poseB_t0, f1i), bodyB.world_vertex(poseB_t0, f2i)));

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    Pose poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = NUM_PIECES; i > 0; i--) {
        ts.push(i / double(NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        Pose poseA_ti1 = Pose::interpolate(poseA_t0, poseA_t1, ti1);
        Pose poseB_ti1 = Pose::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(ipc::point_triangle_distance(
            bodyA.world_vertex(poseA_ti0, vertex_id),
            bodyB.world_vertex(poseB_ti0, f0i),
            bodyB.world_vertex(poseB_ti0, f1i),
            bodyB.world_vertex(poseB_ti0, f2i)));

        if ((distance_ti0 < 0.2 * distance_t0) && ti0 >= 1e-6) {
            toi = ti0;
            is_impacting = true;
            break;
        }

        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval t(ti0, ti1);

        PoseI poseIA_t0 = poseA_t0.cast<Interval>();
        PoseI poseIA_t1 = poseA_t1.cast<Interval>();
        PoseI poseIB_t0 = poseB_t0.cast<Interval>();
        PoseI poseIB_t1 = poseB_t1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_t0, poseIA_t1, t);
        PoseI poseIB = PoseI::interpolate(poseIB_t0, poseIB_t1, t);

        Eigen::MatrixXX3I RA = poseIA.construct_rotation_matrix();
        Eigen::MatrixXX3I RB = poseIB.construct_rotation_matrix();

        using Vector3I = Eigen::Vector3I;

        Vector3I v_t0 = bodyA.world_vertex(poseA_t0, vi).cast<Interval>();
        Vector3I v_t1 = bodyA.world_vertex(poseA_t1, vi).cast<Interval>();
        Vector3I v = bodyA.world_vertex(RA, poseIA.position, vi);
        Interval d = (v - ((v_t1 - v_t0) * t + v_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I f0_t0 = bodyB.world_vertex(poseB_t0, f0i).cast<Interval>();
        Vector3I f0_t1 = bodyB.world_vertex(poseB_t1, f0i).cast<Interval>();
        Vector3I f0 = bodyB.world_vertex(RB, poseIB.position, f0i);
        d = (f0 - ((f0_t1 - f0_t0) * t + f0_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I f1_t0 = bodyB.world_vertex(poseB_t0, f1i).cast<Interval>();
        Vector3I f1_t1 = bodyB.world_vertex(poseB_t1, f1i).cast<Interval>();
        Vector3I f1 = bodyB.world_vertex(RB, poseIB.position, f1i);
        d = (f1 - ((f1_t1 - f1_t0) * t + f1_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        Vector3I f2_t0 = bodyB.world_vertex(poseB_t0, f2i).cast<Interval>();
        Vector3I f2_t1 = bodyB.world_vertex(poseB_t1, f2i).cast<Interval>();
        Vector3I f2 = bodyB.world_vertex(RB, poseIB.position, f2i);
        d = (f2 - ((f2_t1 - f2_t0) * t + f2_t0)).norm();
        min_distance = std::max(min_distance, d.upper());

        if (min_distance >= 0.5 * distance_ti0) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        // spdlog::critical("min_distance={:g}", min_distance);

        double output_tolerance;
        // 0: normal ccd method which only checks t = [0,1]
        // 1: ccd with max_itr and t=[0, t_max]
        const int CCD_TYPE = 1;
        is_impacting = inclusion_ccd::vertexFaceCCD_double(
            bodyA.world_vertex(poseA_ti0, vi),
            bodyB.world_vertex(poseB_ti0, f0i),
            bodyB.world_vertex(poseB_ti0, f1i),
            bodyB.world_vertex(poseB_ti0, f2i),
            bodyA.world_vertex(poseA_ti1, vi),
            bodyB.world_vertex(poseB_ti1, f0i),
            bodyB.world_vertex(poseB_ti1, f1i),
            bodyB.world_vertex(poseB_ti1, f2i),
            { { -1, -1, -1 } }, // rounding error
            min_distance,       // minimum separation distance
            toi,                // time of impact
            1e-6,               // delta
            1.0,                // Maximum time to check
            1e6,                // Maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            break;
        }

        ts.pop();
        ti0 = ti1;
        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }
    // spdlog::trace("vf_ccd_num_subdivision={:d}", num_subdivisions);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    fmt::print("VF {:.16g}\n", timer.getElapsedTime());
#endif

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ccd
