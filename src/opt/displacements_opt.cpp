#include <opt/displacements_opt.hpp>
#include <opt/displacements_opt_nlopt.hpp>

#include <ccd/not_implemented_error.hpp>
#include <ccd/prune_impacts.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    // Optimize the displacments with the volume constraint C(U) ≤ 0.
    double displacements_optimization(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const OptimizationMethod opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt)
    {
        switch (opt_method) {
        case MMA:
        case SLSQP:
            // Both of these methods are in NLopt
            return displacements_optimization_nlopt(V, U, E, volume_epsilon,
                ccd_detection_method,
                opt_method == MMA ? nlopt::LD_MMA : nlopt::LD_SLSQP, max_iter,
                Uopt);
        case IP:
            // Implemented in Ipopt
            throw NotImplementedError(
                "Interior Point optimization method not implemented yet.");
        }
    }

    // Run the entire CCD pipeline.
    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map)
    {
        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(
            V, U, E, ev_impacts, detection_method);
        convert_edge_vertex_to_edge_edge_impacts(E, ev_impacts, ee_impacts);
        ccd::prune_impacts(ee_impacts, edge_impact_map);
    }

}
}
