// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/alternative_displacement_opt.hpp>

#include <autodiff/finitediff.hpp>
#include <ccd/collision_penalty_diff.hpp>
#include <logger.hpp>
#include <opt/barrier_newton_solver.hpp>

namespace ccd {
namespace opt {
    namespace alt {

        // Create a OptimizationProblem for displacment optimization
        void setup_displacement_optimization_problem(const Eigen::MatrixX2d& V,
            const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
            const DetectionMethod ccd_detection_method,
            const bool recompute_collision_set, OptimizationProblem& problem)
        {
            // Detect initial collisions
            // we will keep this set fixed during the optimization
            // -----------------------------------------------------------------
            EdgeEdgeImpacts ee_impacts;
            Eigen::VectorXi edge_impact_map(E.rows());
            detect_collisions(
                V, U, E, ccd_detection_method, ee_impacts, edge_impact_map);
            double starting_epsilon = 2e19;
            for (long i = 0; i < edge_impact_map.size(); i++) {
                starting_epsilon = (edge_impact_map(i) >= 0
                                       && ee_impacts[edge_impact_map(i)].time
                                           < starting_epsilon)
                    ? (ee_impacts[edge_impact_map(i)].time)
                    : (starting_epsilon);
            }
            if (starting_epsilon > 1.0) {
                starting_epsilon = 0.0;
            }

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
            problem.f = [u_](const Eigen::VectorXd& x) {
                return (x - u_).squaredNorm() / 2.0;
            };

            problem.grad_f = [u_](const Eigen::VectorXd& x) { return x - u_; };

            problem.hessian_f = [u_](const Eigen::VectorXd& x) {
                return Eigen::MatrixXd::Identity(x.size(), x.size());
            };

            // g(.) = τ_I - t
            problem.g = [=](const Eigen::VectorXd& x) -> Eigen::VectorXd {
                Eigen::MatrixXd Uk = x;
                Uk.resize(x.rows() / 2, 2);

                Eigen::VectorXd gx;
                if (recompute_collision_set) {
                    EdgeEdgeImpacts new_ee_impacts;
                    Eigen::VectorXi new_edge_impact_map(num_constraints);
                    detect_collisions(V, Uk, E, ccd_detection_method,
                        new_ee_impacts, new_edge_impact_map);
                    ccd::autodiff::compute_penalties_refresh_toi(
                        V, Uk, E, new_ee_impacts, new_edge_impact_map, gx);
                } else {
                    ccd::autodiff::compute_penalties_refresh_toi(
                        V, Uk, E, ee_impacts, edge_impact_map, gx);
                }

                return gx;
            };

            problem.jac_g = [=](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
                Eigen::MatrixXd Uk = x;
                Uk.resize(x.rows() / 2, 2);

                Eigen::MatrixXd jac_gx(num_vars, num_constraints);
                if (recompute_collision_set) {
                    EdgeEdgeImpacts new_ee_impacts;
                    Eigen::VectorXi new_edge_impact_map(num_constraints);
                    detect_collisions(V, Uk, E, ccd_detection_method,
                        new_ee_impacts, new_edge_impact_map);
                    ccd::autodiff::compute_penalties_gradient(
                        V, Uk, E, new_ee_impacts, new_edge_impact_map, jac_gx);
                } else {
                    ccd::autodiff::compute_penalties_gradient(
                        V, Uk, E, ee_impacts, edge_impact_map, jac_gx);
                }
                jac_gx.transposeInPlace();

#ifdef WITH_DERIVATIVE_CHECK
                Eigen::MatrixXd fd_jac_gx;
                ccd::finite_jacobian(x, problem.g, fd_jac_gx);
                if (!ccd::compare_jacobian(jac_gx, fd_jac_gx)) {
                    spdlog::warn("Displ Optimization CHECK_GRADIENT FAILED\n "
                                 "\tgrad={}\n"
                                 "\tfd  ={}",
                        ccd::log::fmt_eigen(jac_gx),
                        ccd::log::fmt_eigen(fd_jac_gx));
                }
#endif
                return jac_gx;
            };

            problem.hessian_g =
                [=](const Eigen::VectorXd& x) -> std::vector<Eigen::MatrixXd> {
                Eigen::MatrixXd Uk = x;
                Uk.resize(x.rows() / 2, 2);

                std::vector<Eigen::MatrixXd> hessian_gx;
                if (recompute_collision_set) {
                    EdgeEdgeImpacts new_ee_impacts;
                    Eigen::VectorXi new_edge_impact_map(num_constraints);
                    detect_collisions(V, Uk, E, ccd_detection_method,
                        new_ee_impacts, new_edge_impact_map);
                    ccd::autodiff::compute_penalties_hessian(V, Uk, E,
                        new_ee_impacts, new_edge_impact_map, hessian_gx);
                } else {
                    ccd::autodiff::compute_penalties_hessian(
                        V, Uk, E, ee_impacts, edge_impact_map, hessian_gx);
                }
                assert(hessian_gx.size() == unsigned(num_constraints));
                assert(
                    hessian_gx.size() == 0 || hessian_gx[0].rows() == num_vars);
                assert(
                    hessian_gx.size() == 0 || hessian_gx[0].cols() == num_vars);
                return hessian_gx;
            };
        }

        // Optimize the displacment opt problem with the given method and
        // starting value.
        OptimizationResults displacement_optimization(
            OptimizationProblem& problem, const Eigen::MatrixX2d& U0,
            SolverSettings& settings)
        {
            // initial value
            Eigen::MatrixXd x0 = U0;
            x0.resize(U0.size(), 1); // Flatten displacements
            problem.x0 = x0;

            OptimizationResults results;

            results = solve_problem_with_barrier_newton(problem, settings);

            results.x.resize(U0.rows(), 2); // Unflatten displacments
            return results;
        }

    } // namespace alt
} // namespace opt
} // namespace ccd
