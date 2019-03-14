// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/solver.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <opt/nlopt_solver.hpp>
#ifdef BUILD_WITH_IPOPT
#include <opt/ipopt_solver.hpp>
#endif
#ifdef BUILD_WITH_OSQP
#include <opt/linearized_constraint_solver.hpp>
#endif
// #include <opt/ncp_solver.hpp> // TODO: Create this file

#include <ccd/not_implemented_error.hpp>
#include <ccd/prune_impacts.hpp>

namespace ccd {
namespace opt {

    // Optimize the displacments with the volume constraint C(U) ≤ 0.
    bool solve_problem(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const double volume_epsilon,
        const DetectionMethod ccd_detection_method,
        const OptimizationMethod opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt)
    {
        switch (opt_method) {
        case MMA:
        case SLSQP:
            // Both of these methods are in NLopt
            return solve_problem_with_nlopt(V, U, E, volume_epsilon,
                ccd_detection_method,
                opt_method == MMA ? nlopt::LD_MMA : nlopt::LD_SLSQP, max_iter,
                Uopt);
        case IP:
#ifdef BUILD_WITH_IPOPT
            // Implemented in Ipopt
            return solve_problem_with_ipopt(
                V, U, E, volume_epsilon, ccd_detection_method, max_iter, Uopt);
#else
            throw NotImplementedError("IPOPT not Enabled");
#endif
        case LINEARIZED_CONSTRAINTS:
#ifdef BUILD_WITH_OSQP
            // Implemented in Ipopt
            return solve_problem_with_linearized_constraints(
                V, U, E, volume_epsilon, ccd_detection_method, max_iter, Uopt);
#else
            throw NotImplementedError("OSQP not Enabled");
#endif
        case NCP:
            throw NotImplementedError(
                "Nonlinear complementarity problem not implemented");
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

    // Save JSON file of optimization objectives per iteration.
    void export_intermediate(const OptimizationMethod method,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints)
    {
        using nlohmann::json;

        std::vector<int> it(objectives.size());
        std::iota(it.begin(), it.end(), 0); // Initalize it with range(0, n)

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
            out_file = std::ofstream("./ip-opt-steps.json");
            std::cout << "./ip-opt-steps.json" << std::endl;
            break;
        case LINEARIZED_CONSTRAINTS:
            out_file = std::ofstream("./linearized-constraints-opt-steps.json");
            std::cout << "./linearized-constraints-opt-steps.json" << std::endl;
            break;
        case NCP:
            out_file = std::ofstream("./ncp-opt-steps.json");
            std::cout << "./ncp-opt-steps.json" << std::endl;
            break;
        }

        out_file << std::setw(4) << figure << std::endl;
    }
} // namespace opt
} // namespace ccd
