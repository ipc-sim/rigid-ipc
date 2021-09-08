// Time-of-impact computation for rigid bodies with angular trajectories.
#include "time_of_impact.hpp"

#include <stack>

#include <tight_inclusion/inclusion_ccd.hpp>

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

// #define TIME_CCD_QUERIES
#ifdef TIME_CCD_QUERIES
#include <igl/Timer.h>
#endif

#include <ipc/ipc.hpp>

#include <ccd/linear/edge_vertex_ccd.hpp>
#include <ccd/rigid/rigid_trajectory_aabb.hpp>
#include <interval/interval.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

// #define USE_FIXED_PIECES
#define USE_DECREASING_DISTANCE_CHECK

namespace ipc::rigid {

typedef Pose<Interval> PoseI;

#ifdef USE_FIXED_PIECES
static const size_t FIXED_NUM_PIECES = 100;
#endif

static const size_t MAX_NUM_SUBDIVISIONS = 1e3;

static const size_t LINEAR_CCD_MAX_ITERATIONS = 1e6;
static const double LINEAR_CCD_TOL = 1e-6;

static const double TRAJECTORY_DISTANCE_FACTOR = 0.5;

#ifdef USE_DECREASING_DISTANCE_CHECK
static const double DECREASING_DISTANCE_FACTOR = 0.2;
static const double DECREASING_DISTANCE_MIN_TIME = 1e-6;
#endif

// 0: normal ccd method which only checks t = [0,1]
// 1: ccd with max_itr and t=[0, t_max]
static const int TIGHT_INCLUSION_CCD_TYPE = 1;

// If true, refine until the toi is not zero. We do not want this because our
// outer CCD algorithmn will perform the refinement in a smarter way.
static const bool TIGHT_INCLUSION_NO_ZERO_TOI = false;

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

bool tight_inclusion_point_edge_ccd(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    const std::array<double, 3>& err,
    const double min_distance,
    double& toi,
    const double tolerance,
    const double tmax,
    const size_t max_iterations,
    double& output_tolerance,
    const int ccd_type = TIGHT_INCLUSION_CCD_TYPE,
    bool no_zero_toi = TIGHT_INCLUSION_NO_ZERO_TOI)
{
    Eigen::Vector3d p_t0_3D, e0_t0_3D, e1_t0_3D, p_t1_3D, e0_t1_3D, e1_t1_3D;
    p_t0_3D.head<2>() = p_t0;
    p_t0_3D.z() = 0;
    e0_t0_3D.head<2>() = e0_t0;
    e0_t0_3D.z() = 0;
    e1_t0_3D.head<2>() = e1_t0;
    e1_t0_3D.z() = 0;
    p_t1_3D.head<2>() = p_t1;
    p_t1_3D.z() = 0;
    e0_t1_3D.head<2>() = e0_t1;
    e0_t1_3D.z() = 0;
    e1_t1_3D.head<2>() = e1_t1;
    e1_t1_3D.z() = 0;

    // NOTE: Use degenerate edge-edge
    return inclusion_ccd::edgeEdgeCCD_double(
        p_t0_3D, p_t0_3D, e0_t0_3D, e1_t0_3D, //
        p_t1_3D, p_t1_3D, e0_t1_3D, e1_t1_3D,
        err,              // rounding error (auto)
        min_distance,     // minimum separation distance
        toi,              // time of impact
        tolerance,        // delta
        tmax,             // maximum time to check
        max_iterations,   // maximum number of iterations
        output_tolerance, // delta_actual
        ccd_type,
        no_zero_toi); // refine until the toi is not zero
}

/// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_vertex_time_of_impact(
    const RigidBody& bodyA, // Body of the vertex
    const PoseD& poseA_t0,  // Pose of bodyA at t=0
    const PoseD& poseA_t1,  // Pose of bodyA at t=1
    size_t vertex_id,       // In bodyA
    const RigidBody& bodyB, // Body of the edge
    const PoseD& poseB_t0,  // Pose of bodyB at t=0
    const PoseD& poseB_t1,  // Pose of bodyB at t=1
    size_t edge_id,         // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double minimum_separation_distance,
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 2);
    assert(minimum_separation_distance >= 0);

    const long vi = vertex_id;
    const long e0i = bodyB.edges(edge_id, 0);
    const long e1i = bodyB.edges(edge_id, 1);

    double distance_t0 = sqrt(point_edge_distance(
        bodyA.world_vertex(poseA_t0, vi), //
        bodyB.world_vertex(poseB_t0, e0i), bodyB.world_vertex(poseB_t0, e1i)));
    if (distance_t0 <= minimum_separation_distance) {
        spdlog::warn(
            "initial distance in edge-vertex CCD is less than MS={:g}!",
            minimum_separation_distance);
        toi = 0;
        return true;
    }

    bool is_impacting = false;
    PoseD poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        PoseD poseA_ti1 = PoseD::interpolate(poseA_t0, poseA_t1, ti1);
        PoseD poseB_ti1 = PoseD::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(point_edge_distance(
            bodyA.world_vertex(poseA_ti1, vi),
            bodyB.world_vertex(poseB_ti1, e0i),
            bodyB.world_vertex(poseB_ti1, e1i)));

#ifdef USE_DECREASING_DISTANCE_CHECK
        if ((distance_ti0 < DECREASING_DISTANCE_FACTOR * distance_t0)
            && ti0 >= DECREASING_DISTANCE_MIN_TIME) {
            toi = ti0;
            is_impacting = true;
            break;
        }
#endif
        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval ti(0, 1);

        PoseI poseIA_ti0 = poseA_ti0.cast<Interval>();
        PoseI poseIA_ti1 = poseA_ti1.cast<Interval>();
        PoseI poseIB_ti0 = poseB_ti0.cast<Interval>();
        PoseI poseIB_ti1 = poseB_ti1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_ti0, poseIA_ti1, ti);
        PoseI poseIB = PoseI::interpolate(poseIB_ti0, poseIB_ti1, ti);

        MatrixMax3I RA = poseIA.construct_rotation_matrix();
        MatrixMax3I RB = poseIB.construct_rotation_matrix();

        Vector2I v_ti0 = bodyA.world_vertex(poseA_ti0, vi).cast<Interval>();
        Vector2I v_ti1 = bodyA.world_vertex(poseA_ti1, vi).cast<Interval>();
        Vector2I v = bodyA.world_vertex(RA, poseIA.position, vi);
        Interval d = (v - ((v_ti1 - v_ti0) * ti + v_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        double v_min_distance = d.upper();

        double e_min_distance = 0;
        Vector2I e0_ti0 = bodyB.world_vertex(poseB_ti0, e0i).cast<Interval>();
        Vector2I e0_ti1 = bodyB.world_vertex(poseB_ti1, e0i).cast<Interval>();
        Vector2I e0 = bodyB.world_vertex(RB, poseIB.position, e0i);
        d = (e0 - ((e0_ti1 - e0_ti0) * ti + e0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        e_min_distance = std::max(e_min_distance, d.upper());

        Vector2I e1_ti0 = bodyB.world_vertex(poseB_ti0, e1i).cast<Interval>();
        Vector2I e1_ti1 = bodyB.world_vertex(poseB_ti1, e1i).cast<Interval>();
        Vector2I e1 = bodyB.world_vertex(RB, poseIB.position, e1i);
        d = (e1 - ((e1_ti1 - e1_ti0) * ti + e1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        e_min_distance = std::max(e_min_distance, d.upper());

        min_distance = v_min_distance + e_min_distance;

        if (min_distance >= TRAJECTORY_DISTANCE_FACTOR * distance_ti0
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        min_distance += minimum_separation_distance;

        double output_tolerance;
        is_impacting = tight_inclusion_point_edge_ccd(
            bodyA.world_vertex(poseA_ti0, vi),
            bodyB.world_vertex(poseB_ti0, e0i),
            bodyB.world_vertex(poseB_ti0, e1i),
            bodyA.world_vertex(poseA_ti1, vi),
            bodyB.world_vertex(poseB_ti1, e0i),
            bodyB.world_vertex(poseB_ti1, e1i),
            { { -1, -1, -1 } },        // rounding error
            min_distance,              // minimum separation distance
            toi,                       // time of impact
            LINEAR_CCD_TOL,            // delta
            1.0,                       // Maximum time to check
            LINEAR_CCD_MAX_ITERATIONS, // Maximum number of iterations
            output_tolerance,          // delta_actual
            TIGHT_INCLUSION_CCD_TYPE,  //
            TIGHT_INCLUSION_NO_ZERO_TOI);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > MS_DIST
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                // spdlog::warn(
                //     "failure=\"Edge-vertex MSCCD says toi=0, but "
                //     "distance_t0={0:g}\" failsafe=\"spliting [{1:g}, "
                //     "{2:g}] into ([{1:g}, {3:g}], [{3:g}, {2:g}])\"",
                //     distance_t0, ti0, ti1, ts.top());
                continue;
            }
            break;
        }

        ts.pop();
        ti0 = ti1;
        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }

    // This time of impact is very dangerous for convergence
    assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_edge_time_of_impact(
    const RigidBody& bodyA, // Body of the first edge
    const PoseD& poseA_t0,  // Pose of bodyA at t=0
    const PoseD& poseA_t1,  // Pose of bodyA at t=1
    size_t edgeA_id,        // In bodyA
    const RigidBody& bodyB, // Body of the second edge
    const PoseD& poseB_t0,  // Pose of bodyB at t=0
    const PoseD& poseB_t1,  // Pose of bodyB at t=1
    size_t edgeB_id,        // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double minimum_separation_distance,
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);
    assert(minimum_separation_distance >= 0);

    const long ea0i = bodyA.edges(edgeA_id, 0);
    const long ea1i = bodyA.edges(edgeA_id, 1);
    const long eb0i = bodyB.edges(edgeB_id, 0);
    const long eb1i = bodyB.edges(edgeB_id, 1);

    double distance_t0 = sqrt(edge_edge_distance(
        bodyA.world_vertex(poseA_t0, ea0i), bodyA.world_vertex(poseA_t0, ea1i),
        bodyB.world_vertex(poseB_t0, eb0i),
        bodyB.world_vertex(poseB_t0, eb1i)));
    if (distance_t0 <= minimum_separation_distance) {
        spdlog::warn(
            "initial distance in edge-edge CCD is less than MS={:g}!",
            minimum_separation_distance);
        toi = 0;
        return true;
    }
    assert(distance_t0 > minimum_separation_distance);

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    PoseD poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        PoseD poseA_ti1 = PoseD::interpolate(poseA_t0, poseA_t1, ti1);
        PoseD poseB_ti1 = PoseD::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(edge_edge_distance(
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 1))));

#ifdef USE_DECREASING_DISTANCE_CHECK
        if (distance_ti0 < DECREASING_DISTANCE_FACTOR * distance_t0
            && ti0 >= DECREASING_DISTANCE_MIN_TIME) {
            toi = ti0;
            is_impacting = true;
            break;
        }
#endif

        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval ti(0, 1);

        PoseI poseIA_ti0 = poseA_ti0.cast<Interval>();
        PoseI poseIA_ti1 = poseA_ti1.cast<Interval>();
        PoseI poseIB_ti0 = poseB_ti0.cast<Interval>();
        PoseI poseIB_ti1 = poseB_ti1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_ti0, poseIA_ti1, ti);
        PoseI poseIB = PoseI::interpolate(poseIB_ti0, poseIB_ti1, ti);

        MatrixMax3I RA = poseIA.construct_rotation_matrix();
        MatrixMax3I RB = poseIB.construct_rotation_matrix();

        double min_ea_distance = 0;
        Vector3I ea0_ti0 = bodyA.world_vertex(poseA_ti0, ea0i).cast<Interval>();
        Vector3I ea0_ti1 = bodyA.world_vertex(poseA_ti1, ea0i).cast<Interval>();
        Vector3I ea0 = bodyA.world_vertex(RA, poseIA.position, ea0i);
        Interval d = (ea0 - ((ea0_ti1 - ea0_ti0) * ti + ea0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_ea_distance = std::max(min_ea_distance, d.upper());

        Vector3I ea1_ti0 = bodyA.world_vertex(poseA_ti0, ea1i).cast<Interval>();
        Vector3I ea1_ti1 = bodyA.world_vertex(poseA_ti1, ea1i).cast<Interval>();
        Vector3I ea1 = bodyA.world_vertex(RA, poseIA.position, ea1i);
        d = (ea1 - ((ea1_ti1 - ea1_ti0) * ti + ea1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_ea_distance = std::max(min_ea_distance, d.upper());

        double min_eb_distance = 0;
        Vector3I eb0_ti0 = bodyB.world_vertex(poseB_ti0, eb0i).cast<Interval>();
        Vector3I eb0_ti1 = bodyB.world_vertex(poseB_ti1, eb0i).cast<Interval>();
        Vector3I eb0 = bodyB.world_vertex(RB, poseIB.position, eb0i);
        d = (eb0 - ((eb0_ti1 - eb0_ti0) * ti + eb0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_eb_distance = std::max(min_eb_distance, d.upper());

        Vector3I eb1_ti0 = bodyB.world_vertex(poseB_ti0, eb1i).cast<Interval>();
        Vector3I eb1_ti1 = bodyB.world_vertex(poseB_ti1, eb1i).cast<Interval>();
        Vector3I eb1 = bodyB.world_vertex(RB, poseIB.position, eb1i);
        d = (eb1 - ((eb1_ti1 - eb1_ti0) * ti + eb1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_eb_distance = std::max(min_eb_distance, d.upper());

        min_distance = min_ea_distance + min_eb_distance;

        if (min_distance >= TRAJECTORY_DISTANCE_FACTOR * distance_ti0
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        min_distance += minimum_separation_distance;
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
            { { -1, -1, -1 } },        // rounding error
            min_distance,              // minimum separation distance
            toi,                       // time of impact
            LINEAR_CCD_TOL,            // delta
            1.0,                       // Maximum time to check
            LINEAR_CCD_MAX_ITERATIONS, // Maximum number of iterations
            output_tolerance,          // delta_actual
            TIGHT_INCLUSION_CCD_TYPE,  //
            TIGHT_INCLUSION_NO_ZERO_TOI);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > MS_DIST
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                // spdlog::warn(
                //     "failure=\"Edge-edge MSCCD says toi=0, but "
                //     "distance_t0={0:g}\" failsafe=\"spliting [{1:g}, "
                //     "{2:g}] into ([{1:g}, {3:g}], [{3:g}, {2:g}])\"",
                //     distance_t0, ti0, ti1, ts.top());
                continue;
            }
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
    assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Face-Vertex

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_face_vertex_time_of_impact(
    const RigidBody& bodyA, // Body of the vertex
    const PoseD& poseA_t0,  // Pose of bodyA at t=0
    const PoseD& poseA_t1,  // Pose of bodyA at t=1
    size_t vertex_id,       // In bodyA
    const RigidBody& bodyB, // Body of the triangle
    const PoseD& poseB_t0,  // Pose of bodyB at t=0
    const PoseD& poseB_t1,  // Pose of bodyB at t=1
    size_t face_id,         // In bodyB
    double& toi,
    double earliest_toi, // Only search for collision in [0, earliest_toi]
    double minimum_separation_distance,
    double toi_tolerance)
{
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);
    assert(minimum_separation_distance >= 0);

    const long vi = vertex_id;
    const long f0i = bodyB.faces(face_id, 0);
    const long f1i = bodyB.faces(face_id, 1);
    const long f2i = bodyB.faces(face_id, 2);

    double distance_t0 = sqrt(point_triangle_distance(
        bodyA.world_vertex(poseA_t0, vi), bodyB.world_vertex(poseB_t0, f0i),
        bodyB.world_vertex(poseB_t0, f1i), bodyB.world_vertex(poseB_t0, f2i)));
    if (distance_t0 <= minimum_separation_distance) {
        spdlog::warn(
            "initial distance in faces-vertex CCD is less than MS={:g}!",
            minimum_separation_distance);
        toi = 0;
        return true;
    }

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    PoseD poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        PoseD poseA_ti1 = PoseD::interpolate(poseA_t0, poseA_t1, ti1);
        PoseD poseB_ti1 = PoseD::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(point_triangle_distance(
            bodyA.world_vertex(poseA_ti0, vi),
            bodyB.world_vertex(poseB_ti0, f0i),
            bodyB.world_vertex(poseB_ti0, f1i),
            bodyB.world_vertex(poseB_ti0, f2i)));

#ifdef USE_DECREASING_DISTANCE_CHECK
        if ((distance_ti0 < DECREASING_DISTANCE_FACTOR * distance_t0)
            && ti0 >= DECREASING_DISTANCE_MIN_TIME) {
            toi = ti0;
            is_impacting = true;
            break;
        }
#endif
        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval ti(0, 1);

        PoseI poseIA_ti0 = poseA_ti0.cast<Interval>();
        PoseI poseIA_ti1 = poseA_ti1.cast<Interval>();
        PoseI poseIB_ti0 = poseB_ti0.cast<Interval>();
        PoseI poseIB_ti1 = poseB_ti1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_ti0, poseIA_ti1, ti);
        PoseI poseIB = PoseI::interpolate(poseIB_ti0, poseIB_ti1, ti);

        MatrixMax3I RA = poseIA.construct_rotation_matrix();
        MatrixMax3I RB = poseIB.construct_rotation_matrix();

        Vector3I v_ti0 = bodyA.world_vertex(poseA_ti0, vi).cast<Interval>();
        Vector3I v_ti1 = bodyA.world_vertex(poseA_ti1, vi).cast<Interval>();
        Vector3I v = bodyA.world_vertex(RA, poseIA.position, vi);
        Interval d = (v - ((v_ti1 - v_ti0) * ti + v_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        double v_min_distance = d.upper();

        double f_min_distance = 0;
        Vector3I f0_ti0 = bodyB.world_vertex(poseB_ti0, f0i).cast<Interval>();
        Vector3I f0_ti1 = bodyB.world_vertex(poseB_ti1, f0i).cast<Interval>();
        Vector3I f0 = bodyB.world_vertex(RB, poseIB.position, f0i);
        d = (f0 - ((f0_ti1 - f0_ti0) * ti + f0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        f_min_distance = std::max(f_min_distance, d.upper());

        Vector3I f1_ti0 = bodyB.world_vertex(poseB_ti0, f1i).cast<Interval>();
        Vector3I f1_ti1 = bodyB.world_vertex(poseB_ti1, f1i).cast<Interval>();
        Vector3I f1 = bodyB.world_vertex(RB, poseIB.position, f1i);
        d = (f1 - ((f1_ti1 - f1_ti0) * ti + f1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        f_min_distance = std::max(f_min_distance, d.upper());

        Vector3I f2_ti0 = bodyB.world_vertex(poseB_ti0, f2i).cast<Interval>();
        Vector3I f2_ti1 = bodyB.world_vertex(poseB_ti1, f2i).cast<Interval>();
        Vector3I f2 = bodyB.world_vertex(RB, poseIB.position, f2i);
        d = (f2 - ((f2_ti1 - f2_ti0) * ti + f2_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        f_min_distance = std::max(f_min_distance, d.upper());

        min_distance = v_min_distance + f_min_distance;

        if (min_distance >= TRAJECTORY_DISTANCE_FACTOR * distance_ti0
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        min_distance += minimum_separation_distance;

        double output_tolerance;
        is_impacting = inclusion_ccd::vertexFaceCCD_double(
            bodyA.world_vertex(poseA_ti0, vi),
            bodyB.world_vertex(poseB_ti0, f0i),
            bodyB.world_vertex(poseB_ti0, f1i),
            bodyB.world_vertex(poseB_ti0, f2i),
            bodyA.world_vertex(poseA_ti1, vi),
            bodyB.world_vertex(poseB_ti1, f0i),
            bodyB.world_vertex(poseB_ti1, f1i),
            bodyB.world_vertex(poseB_ti1, f2i),
            { { -1, -1, -1 } },        // rounding error
            min_distance,              // minimum separation distance
            toi,                       // time of impact
            LINEAR_CCD_TOL,            // delta
            1.0,                       // Maximum time to check
            LINEAR_CCD_MAX_ITERATIONS, // Maximum number of iterations
            output_tolerance,          // delta_actual
            TIGHT_INCLUSION_CCD_TYPE,  //
            TIGHT_INCLUSION_NO_ZERO_TOI);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > MS_DIST
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                // spdlog::warn(
                //     "failure=\"Face-vertex MSCCD says toi=0, but "
                //     "distance_t0={0:g}\" failsafe=\"spliting [{1:g}, "
                //     "{2:g}] into ([{1:g}, {3:g}], [{3:g}, {2:g}])\"",
                //     distance_t0, ti0, ti1, ts.top());
                continue;
            }
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
    assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ipc::rigid
