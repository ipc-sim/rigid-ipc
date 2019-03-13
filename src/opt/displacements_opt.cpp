#include <opt/displacements_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <opt/displacements_opt_nlopt.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/displacements_opt_ipopt.hpp>
#endif

#include <ccd/not_implemented_error.hpp>
#include <ccd/prune_impacts.hpp>

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
#ifdef BUILD_WITH_IPOPT
        case IP:
            // Implemented in Ipopt
            return displacements_optimization_ipopt(
                V, U, E, volume_epsilon, ccd_detection_method, max_iter, Uopt);
#else
        default:
            throw NotImplementedError("IPOPT not Enabled");
#endif
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

    void export_intermediate(const OptimizationMethod method,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints)
    {
        using nlohmann::json;

        std::vector<int> it(objectives.size());
        std::iota(it.begin(), it.end(), 0);

        json data;
        data["x"] = it;
        data["objectives"] = objectives;
        data["constraints"] = constraints;
        json figure;
        figure["data"] = data;
        std::ofstream out_file;
        switch (method) {
        case MMA:
            out_file = std::ofstream(
                "../figures/optimization-iterations/mma-opt-steps.json");
            break;
        case SLSQP:
            out_file = std::ofstream(
                "../figures/optimization-iterations/mma-opt-steps.json");
            break;
        case IP:
            out_file = std::ofstream(
                "./ip-opt-steps.json");
             std::cout << "./ip-opt-steps.json" << std::endl;
            break;
        }

        out_file << std::setw(4) << figure << std::endl;
    }
}
}
