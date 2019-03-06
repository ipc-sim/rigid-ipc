#include <opt/displacements_nlopt.hpp>
#include <opt/displacements_opt.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    /**
     * @brief Class for storing optimization data.
     * This is passed to the volume constraints to compute the STIVs.
     */
    class OptData {
    public:
        const Eigen::MatrixX2d& V;   ///< @brief Vertices.
        const Eigen::MatrixX2d& U;   ///< @brief Displacments.
        const Eigen::MatrixX2i& E;   ///< @brief Edges.
        const double volume_epsilon; ///< @brief Epsilon for STIV computation.
        /// @brief Collision detection method
        const DetectionMethod ccd_detection_method;

        /**
         * @brief Store the data needed for optimization.
         * Stores everything by reference because all values should be valid
         * during optimization.
         *
         * @param[in] V Vertices
         * @param[in] U Displacments
         * @param[in] E Edges
         * @param[in] volume_epsilon Epsilon for STIV computation.
         * @param[in] ccd_detection_method Collision detection method.
         */
        OptData(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
            const Eigen::MatrixX2i& E, const double volume_epsilon,
            const DetectionMethod ccd_detection_method)
            : V(V)
            , U(U)
            , E(E)
            , volume_epsilon(volume_epsilon)
            , ccd_detection_method(ccd_detection_method)
        {
        }
    };

    /**
     * @brief Computes the optimization objective (\f$\|U - U_0\|^2\f$).
     *
     * @param[in] u Optimization variable for which to compute \f$f(u)\f$ and
     *              \f$\nabla f(u)\f$.
     * @param[out] grad Location to store \f$\nabla f(x)\f$.
     * @param[in] u_flat A pointer to the flattend initial displacments.
     * @return Returns the value \f$\|U - U_0\|^2\f$.
     */
    double objective(
        const std::vector<double>& u, std::vector<double>& grad, void* u_flat)
    {
        assert(u_flat); // Need the inital displacments to compute anything.
        const Eigen::MatrixXd U_FLAT
            = *static_cast<const Eigen::MatrixXd*>(u_flat);
        // Map the input displacment ot a Eigen Vector.
        const Eigen::VectorXd U
            = Eigen::Map<const Eigen::VectorXd>(u.data(), long(u.size()));
        // Compute the difference between the original and input displacments.
        const Eigen::VectorXd DIFF_U = U - U_FLAT;
        // Only compute the gradient as needed.
        if (!grad.empty()) {
            // This should always be true according to NLopt.
            assert(size_t(grad.size()) == size_t(DIFF_U.size()));
            // ∇f(U) = ∇ (U - U0)^T * (U - U0) = 2 * (U - U0)
            Eigen::VectorXd::Map(&grad[0], long(grad.size())) = 2 * DIFF_U;
        }
        return DIFF_U.squaredNorm(); // f(U) = (U - U0)^T * (U - U0)
    }

    /**
     * @brief Compute the STIV (and ∇ STIV if needed).
     * These are NLopt inequality constraints where volume_constraints(U) ≤ 0.
     *
     * @param[in] m Number of constraints.
     * @param[out] results Array to store the constraint results in.
     * @param[in] n Number of dimensions for the optimization varaible
     *          (n = 2 * |V|).
     * @param[in] u The optimization displacment values.
     * @param[out] grad Array for where to store the volume gradients.
     * @param[in] d
     */
    void volume_constraints(unsigned m, double* results, unsigned n,
        const double* u, double* grad, void* data)
    {
        assert(data);
        OptData* _ = static_cast<OptData*>(data);
        const Eigen::MatrixX2d& V = _->V;
        const Eigen::MatrixX2i& E = _->E;
        const double VOLUME_EPSILON = _->volume_epsilon;
        const DetectionMethod CCD_DETECTION_METHOD = _->ccd_detection_method;

        // Reshape the input displacments to a |V| x 2 matrix.
        Eigen::MatrixXd U = Eigen::Map<const Eigen::VectorXd>(u, n);
        U.resize(U.rows() / 2, 2);

        // Detect the collisions for the input displacments.
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, CCD_DETECTION_METHOD, ee_impacts, edge_impact_map);

        // Compute the STIVs
        Eigen::VectorXd volumes;
        ccd::compute_volumes(
            V, U, E, ee_impacts, edge_impact_map, VOLUME_EPSILON, volumes);

        // Only compute the gradient as needed.
        if (grad) {
            // Gradient of the volumes (n x m)
            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(V, U, E, ee_impacts,
                edge_impact_map, VOLUME_EPSILON, volume_gradient);
            // This should always be true according to NLopt
            assert(size_t(m * n) == size_t(volume_gradient.size()));
            volume_gradient.resize(m * n, 1);
            Eigen::VectorXd::Map(grad, m * n) = -volume_gradient;
        }

        // We want volumes[i] ≥ 0, but NLopt expects constraints(x) ≤ 0,
        // so negate the volume and volume gradient.
        assert(size_t(volumes.size()) == size_t(m));
        // Log the resulting volumes.
        std::cout << "Volumes: " << std::endl;
        for (unsigned i = 0; i < m; i++) {
            std::cout << -volumes[i] << std::endl;
        }
        // Store the volumes in the results.
        Eigen::VectorXd::Map(results, m) = -volumes;
    }

    // Optimize the displacments using NLopt.
    void displacements_nlopt_step(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        Eigen::MatrixX2d& Uopt, const nlopt::algorithm opt_method)
    {
        Eigen::MatrixXd U_flat = U;
        U_flat.resize(U.size(), 1);

        // Construct the data class to pass to the optimization.
        OptData data(V, U, E, volume_epsilon, ccd_detection_method);

        nlopt::opt opt(opt_method, unsigned(U.size())); // Optimization object
        opt.set_min_objective(objective, &U_flat);
        std::vector<double> tol(unsigned(E.rows()), 1e-8); // Tolerances
        // m volume constraints
        opt.add_inequality_mconstraint(volume_constraints, &data, tol);
        opt.set_xtol_rel(1e-12);
        opt.set_maxeval(1000);
        std::vector<double> x(size_t(U.size()), 0); // initial guess is zero
        double minf; // the minimum objective value, upon return
        opt.optimize(x, minf);
        Eigen::MatrixXd Xopt
            = Eigen::Map<Eigen::VectorXd>(x.data(), long(x.size()));
        Xopt.resize(Xopt.rows() / 2, 2);
        Uopt = Xopt;
        std::cout << "Optimal Displacments:\n" << Uopt << std::endl;
    }

}
}
