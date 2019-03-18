#include <opt/displacements_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>
#include <ccd/prune_impacts.hpp>

#include <opt/displacements_opt_nlopt.hpp>
#include <opt/minimize.hpp>

namespace ccd {
namespace opt {

    void setup_displacement_optimization(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        OptimizationProblem& problem)
    {
        // Detect initial collisions
        // we will keep this set fixed during the optimization
        // -----------------------------------------------------------------
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, ccd_detection_method, ee_impacts, edge_impact_map);

        /// Min ||U - Uk||^2
        /// s.t V(U) >= 0
        Eigen::MatrixXd u_ = U;
        u_.resize(U.size(), 1);

        // Setup problem
        // -----------------------------------------------------------------
        int num_vars, num_constraints;
        problem.num_vars = num_vars = int(U.size());
        problem.num_constraints = num_constraints = int(E.rows());

        problem.x0.resize(num_vars);
        problem.x_lower.resize(num_vars);
        problem.x_upper.resize(num_vars);
        problem.g_lower.resize(num_constraints);
        problem.g_upper.resize(num_constraints);

        problem.x_lower.setConstant(NO_LOWER_BOUND);
        problem.x_upper.setConstant(NO_UPPER_BOUND);
        problem.g_lower.setConstant(0.0);
        problem.g_upper.setConstant(NO_UPPER_BOUND);

        // we use lambda by copy [=] to fix the entried to current values

        problem.f = [u_](const Eigen::VectorXd& x) -> double {
            return (x - u_).squaredNorm();
        };

        problem.grad_f = [u_](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return 2.0 * (x - u_);
        };

        problem.g = [=](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);

            Eigen::VectorXd volumes;
            ccd::autodiff::compute_volumes_refresh_toi(
                V, Uk, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);
            return volumes;
        };

        problem.jac_g = [=](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);

            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(V, Uk, E, ee_impacts,
                edge_impact_map, volume_epsilon, volume_gradient);
            assert(volume_gradient.rows() == num_vars);
            assert(volume_gradient.cols() == num_constraints);

            // Note: standard is to be num_constraints x num_vars
            return volume_gradient.transpose();
        };
    }

    OptimizationResult displacements_optimization(
        const OptimizationMethod& opt_method, const Eigen::MatrixX2d& U0,
        OptimizationProblem& problem, std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history, std::vector<double>& g_history)
    {

        problem.intermediate_cb
            = [&u_history, &f_history, &g_history, &problem](
                  const Eigen::VectorXd& x, const double obj_value,
                  const Eigen::VectorXd& /*dual*/,
                  const int /*iteration*/) -> void {
            Eigen::MatrixXd u = x;
            u.resize(x.rows() / 2, 2);

            u_history.push_back(u);
            f_history.push_back(obj_value);
            g_history.push_back(problem.g(x).sum());
        };

        // initial value
        Eigen::MatrixXd x0 = U0;
        x0.resize(U0.size(), 1);
        problem.x0 = x0;

        OptimizationResult result = minimize(opt_method, problem);
        result.x.resize(U0.rows(), 2);

        return result;
    }

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

        default:
            throw NotImplementedError("IPOPT not Enabled");
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

    void export_intermediate(const std::string filename,
        const std::vector<Eigen::Matrix2d>& displacements,
        const std::vector<double>& objectives,
        const std::vector<double>& constraints)
    {
        using nlohmann::json;

        std::vector<int> it(objectives.size());
        std::iota(it.begin(), it.end(), 0);

        json figure, data, disp;

        data["x"] = it;
        data["objectives"] = objectives;
        data["constraints"] = constraints;

        figure["data"] = data;

        std::ofstream out_file(filename);
        out_file << std::setw(4) << figure << std::endl;
    }
} // namespace opt
} // namespace ccd
