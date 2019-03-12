#ifdef BUILD_WITH_IPOPT
#include "displacements_opt_ipopt.hpp"

#include <iostream>
#include <nlohmann/json.hpp>

#include <opt/IpEigenInterfaceTNLP.hpp>
#include <opt/displacements_opt.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

namespace ccd {
namespace opt {
    double displacements_optimization_ipopt(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const unsigned max_iter, Eigen::MatrixX2d& Uopt)
    {
        using namespace ccd::opt::ip;
        // Detect initial collisions
        // we will keep this set fixed during the optimization
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, ccd_detection_method, ee_impacts, edge_impact_map);

        /// Min ||U - Uk||^2
        /// s.t V(U) >= 0
        Eigen::MatrixXd u_ = U;
        u_.resize(U.size(), 1);

        int num_vars = int(U.size()), num_constraints = int(E.rows());
        Eigen::VectorXd x_lower(num_vars), x_upper(num_vars);
        Eigen::VectorXd g_lower(num_constraints), g_upper(num_constraints);

        x_lower.setConstant(NO_LOWER_BOUND);
        x_upper.setConstant(NO_UPPER_BOUND);
        g_lower.setConstant(0.0);
        g_upper.setConstant(NO_UPPER_BOUND);

        callback_f f = [u_](const Eigen::VectorXd& x) -> double {
            return (x - u_).squaredNorm();
        };

        callback_grad_f grad_f
            = [u_](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return 2.0 * (x - u_);
        };

        callback_g g = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);

            Eigen::VectorXd volumes;
            ccd::autodiff::compute_volumes_refresh_toi(
                V, Uk, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);
            return volumes;
        };

        callback_jac_g jac_g;
        jac_g = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);

            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(V, Uk, E, ee_impacts,
                edge_impact_map, volume_epsilon, volume_gradient);
            assert(volume_gradient.rows() == num_vars);
            assert(volume_gradient.cols() == num_constraints);

            // the derivative of constraint g^{(i)} with respect to variable
            // x^{(j)} is placed in row i and column j.
            return volume_gradient.transpose();
        };


        std::vector<double> f_history;
        std::vector<double> g_history;
        callback_intermediate callback;
        callback = [&f_history, &g_history, g](const Eigen::VectorXd& x,
                       const double obj_value, const Eigen::VectorXd& /*dual*/,
                       const int iteration) -> void {
            f_history.push_back(obj_value);
            g_history.push_back(g(x).sum());
            std::cout << "len "<< f_history.size() << std::endl;
            std::cout << "it " << iteration << std::endl;
        };

        // initial value
        Eigen::MatrixXd x0 = Uopt;
        x0.resize(U.size(), 1);

        auto result = ccd::opt::ip::minimize_ipopt(f, x0, grad_f, x_lower,
            x_upper, num_constraints, g, jac_g, g_lower, g_upper,
            /*verbosity=*/5, max_iter, callback);

        export_intermediate(IP, f_history, g_history);

        x0 = result.x;
        x0.resize(U.rows(), 2);
        Uopt = x0;
        return result.value;
    }

}

}

#endif
