// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/displacement_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>
#include <ccd/prune_impacts.hpp>

#include <opt/ncp_solver.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

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

    // Create a OptimizationProblem for displacment optimization
    void setup_displacement_optimization_problem(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const bool recompute_collision_set, OptimizationProblem& problem)
    {
        // Detect initial collisions
        // we will keep this set fixed during the optimization
        // -----------------------------------------------------------------
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, ccd_detection_method, ee_impacts, edge_impact_map);

        /// Min 1/2||U - Uk||^2
        /// s.t V(U) >= 0
        Eigen::MatrixXd u_ = U; // U_flat
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
            return (x - u_).squaredNorm() / 2.0;
        };

        problem.grad_f = [u_](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return (x - u_);
        };

        auto g
            = [V, E, volume_epsilon](const Eigen::MatrixXd& Uk,
                  const EdgeEdgeImpacts& ee_impacts,
                  const Eigen::VectorXi& edge_impact_map) -> Eigen::VectorXd {
            Eigen::VectorXd volumes;
            ccd::autodiff::compute_volumes_refresh_toi(
                V, Uk, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);
            return volumes;
        };

        auto jac_g
            = [V, E, volume_epsilon, num_vars, num_constraints](
                  const Eigen::MatrixXd& Uk, const EdgeEdgeImpacts& ee_impacts,
                  const Eigen::VectorXi& edge_impact_map) -> Eigen::MatrixXd {
            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(V, Uk, E, ee_impacts,
                edge_impact_map, volume_epsilon, volume_gradient);
            assert(volume_gradient.rows() == num_vars);
            assert(volume_gradient.cols() == num_constraints);

            // Note: standard is to be num_constraints x num_vars
            return volume_gradient.transpose();
        };
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
        problem.g = [=](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);

            if (recompute_collision_set) {
                EdgeEdgeImpacts new_ee_impacts;
                Eigen::VectorXi new_edge_impact_map(num_constraints);
                detect_collisions(V, Uk, E, ccd_detection_method,
                    new_ee_impacts, new_edge_impact_map);
                return g(Uk, new_ee_impacts, new_edge_impact_map);
            }
            return g(Uk, ee_impacts, edge_impact_map);
        };

        problem.jac_g = [=](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);
            if (recompute_collision_set) {
                EdgeEdgeImpacts new_ee_impacts;
                Eigen::VectorXi new_edge_impact_map(num_constraints);
                detect_collisions(V, Uk, E, ccd_detection_method,
                    new_ee_impacts, new_edge_impact_map);
                return jac_g(Uk, new_ee_impacts, new_edge_impact_map);
            }
            return jac_g(Uk, ee_impacts, edge_impact_map);
        };
#pragma clang diagnostic pop
    }

    // Optimize the displacment opt problem with the given method and starting
    // value.
    OptimizationResults displacement_optimization(OptimizationProblem& problem,
        const Eigen::MatrixX2d& U0, SolverSettings& settings)
    {

        // initial value
        Eigen::MatrixXd x0 = U0;
        x0.resize(U0.size(), 1); // Flatten displacements
        problem.x0 = x0;

        OptimizationResults result;
        if (settings.method == NCP) {
            result = solve_ncp_displacement_optimization(problem, settings);
        } else {
            result = solve_problem(problem, settings);
        }
        result.x.resize(U0.rows(), 2); // Unflatten displacments
        return result;
    } // namespace opt

    OptimizationResults solve_ncp_displacement_optimization(
        OptimizationProblem& problem, SolverSettings& settings)
    {
        // Solves the KKT conditions of the Optimization Problem
        //  (U - Uk) = \nabla g(U)
        //  s.t V(U) >= 0
        int num_vars = int(problem.x0.rows());
        Eigen::SparseMatrix<double> A(num_vars, num_vars);
        A.setIdentity();

        Eigen::VectorXd b, x_opt, lambda_opt;
        // obtain Uk from grad_f using x=0
        b = -problem.grad_f(Eigen::VectorXd::Zero(num_vars));

        int num_it = 0;
        callback_intermediate_ncp callback =
            [&problem, &settings, &num_it](const Eigen::VectorXd& x,
                const Eigen::VectorXd& alpha, const double gamma) {
                settings.intermediate_cb(x, problem.f(x), alpha, gamma, num_it);
                num_it += 1;
            };
        OptimizationResults result;
        result.success = solve_ncp(A, b, problem.g, problem.jac_g,
            settings.max_iter, callback, GRADIENT_ONLY, settings.lcp_solver,
            x_opt, lambda_opt, /*check_convergence=*/false);
        result.x = x_opt;
        result.minf = problem.f(x_opt);

        return result;
    }

} // namespace opt
} // namespace ccd
