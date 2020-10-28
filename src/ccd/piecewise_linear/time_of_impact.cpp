// Time-of-impact computation for rigid bodies with angular trajectories.
#include "time_of_impact.hpp"

#include <tight_inclusion/inclusion_ccd.hpp>

// #define TIME_CCD_QUERIES
#ifdef TIME_CCD_QUERIES
#include <igl/Timer.h>
#endif

#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

static const int NUM_PIECES = 100;

////////////////////////////////////////////////////////////////////////////////
// Edge-Vertex

/// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_vertex_time_of_impact(
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
    throw NotImplementedError(
        "Piecewise linear edge vertex CCD is not implemented yet!");
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_edge_edge_time_of_impact(
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
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);
    assert(dim == 3);

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    physics::Pose<double> poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    for (int i = 1; i <= NUM_PIECES; i++) {
        double ti0 = (i - 1) / double(NUM_PIECES) * earliest_toi;
        double ti1 = i / double(NUM_PIECES) * earliest_toi;
        physics::Pose<double> poseA_ti1 =
            physics::Pose<double>::interpolate(poseA_t0, poseA_t1, ti1);
        physics::Pose<double> poseB_ti1 =
            physics::Pose<double>::interpolate(poseB_t0, poseB_t1, ti1);

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
            0,                  // minimum separation distance
            toi,                // time of impact
            1e-4,               // delta
            1.0,                // Maximum time to check
            -1,                 // Maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            break;
        }

        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }

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

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Face-Vertex

// Find time-of-impact between two rigid bodies
bool compute_piecewise_linear_face_vertex_time_of_impact(
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

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    physics::Pose<double> poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    for (int i = 1; i <= NUM_PIECES; i++) {
        double ti1 = i / double(NUM_PIECES) * earliest_toi;
        physics::Pose<double> poseA_ti1 =
            physics::Pose<double>::interpolate(poseA_t0, poseA_t1, ti1);
        physics::Pose<double> poseB_ti1 =
            physics::Pose<double>::interpolate(poseB_t0, poseB_t1, ti1);

        double output_tolerance;
        // 0: normal ccd method which only checks t = [0,1]
        // 1: ccd with max_itr and t=[0, t_max]
        const int CCD_TYPE = 1;
        is_impacting = inclusion_ccd::vertexFaceCCD_double(
            bodyA.world_vertex(poseA_ti0, vertex_id),
            bodyB.world_vertex(poseB_ti0, bodyB.faces(face_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.faces(face_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.faces(face_id, 2)),
            bodyA.world_vertex(poseA_ti1, vertex_id),
            bodyB.world_vertex(poseB_ti1, bodyB.faces(face_id, 0)),
            bodyB.world_vertex(poseB_ti1, bodyB.faces(face_id, 1)),
            bodyB.world_vertex(poseB_ti1, bodyB.faces(face_id, 2)),
            { { -1, -1, -1 } }, // rounding error
            0,                  // minimum separation distance
            toi,                // time of impact
            1e-6,               // delta
            1.0,                // Maximum time to check
            1e6,                // Maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);

        if (is_impacting) {
            double ti0 = (i - 1) / double(NUM_PIECES) * earliest_toi;
            toi = (ti1 - ti0) * toi + ti0;
            break;
        }

        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }

#ifdef TIME_CCD_QUERIES
    timer.stop();
    std::cout << "VF " << timer.getElapsedTime() << std::endl;
#endif

    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ccd
