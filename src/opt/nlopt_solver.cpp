// Solve a optimization problem with NLopt.
#include <opt/nlopt_solver.hpp>
#include <opt/solver.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    /**
     * @brief Class for storing objective data.
     * This is passed to the objective
     */
    class ObjectiveData {
    public:
        /// @brief Flattened displacment vector, \f$U_0\f$.
        const Eigen::MatrixXd& U_FLAT;
        /// @breif Vector of objectives for every iteration.
        std::vector<double>& objective_per_iteration;

        /**
         * @brief Store the data needed for computing the objective.
         * Stores everything by reference because all values should be valid
         * during optimization.
         *
         * @param[in] U_FLAT Flattened displacments.
         * @param[in] objective_per_iteration Vector of objectives for every
         * iteration.
         */
        ObjectiveData(const Eigen::MatrixXd& U_FLAT,
            std::vector<double>& objective_per_iteration)
            : U_FLAT(U_FLAT)
            , objective_per_iteration(objective_per_iteration)
        {
        }
    };

    /**
     * @brief Class for storing constraint data.
     * This is passed to the volume constraints to compute the STIVs.
     */
    class ConstraintData {
    public:
        const Eigen::MatrixX2d& V;   ///< @brief Vertices.
        const Eigen::MatrixX2i& E;   ///< @brief Edges.
        const double VOLUME_EPSILON; ///< @brief Epsilon for STIV computation.
        /// @brief Collision detection method
        const DetectionMethod CCD_DETECTION_METHOD;
        /// @brief Initial edge-edge impacts.
        const EdgeEdgeImpacts& EE_IMPACTS;
        /// @brief Earliest inital impact per edge.
        const Eigen::VectorXi& EDGE_IMPACT_MAP;
        /// @breif Vector of constraints for every iteration.
        std::vector<double>& sum_constraints_per_iteration;

        /**
         * @brief Store the data needed for optimization.
         * Stores everything by reference because all values should be valid
         * during optimization.
         *
         * @param[in] V Vertices
         * @param[in] U Displacments
         * @param[in] E Edges
         * @param[in] VOLUME_EPSILON Epsilon for STIV computation.
         * @param[in] CCD_DETECTION_METHOD Collision detection method.
         * @param[in] EE_IMPACTS Initial edge-edge impacts.
         * @param[in] EDGE_IMPACT_MAP Earliest inital impact per edge.
         * @param[in] sum_constraints_per_iteration Vector of constraints for
         * every iteration.
         */
        ConstraintData(const Eigen::MatrixX2d& V, const Eigen::MatrixX2i& E,
            const double VOLUME_EPSILON,
            const DetectionMethod CCD_DETECTION_METHOD,
            const EdgeEdgeImpacts& EE_IMPACTS,
            const Eigen::VectorXi& EDGE_IMPACT_MAP,
            std::vector<double>& sum_constraints_per_iteration)
            : V(V)
            , E(E)
            , VOLUME_EPSILON(VOLUME_EPSILON)
            , CCD_DETECTION_METHOD(CCD_DETECTION_METHOD)
            , EE_IMPACTS(EE_IMPACTS)
            , EDGE_IMPACT_MAP(EDGE_IMPACT_MAP)
            , sum_constraints_per_iteration(sum_constraints_per_iteration)
        {
        }
    };

    /**
     * @brief Computes the optimization objective (\f$\|U - U_0\|^2\f$).
     *
     * @param[in] u Optimization variable for which to compute \f$f(u)\f$ and
     *              \f$\nabla f(u)\f$.
     * @param[out] grad Location to store \f$\nabla f(x)\f$.
     * @param[in,out] data A pointer to the objective data.
     * @return Returns the value \f$\|U - U_0\|^2\f$.
     */
    double objective(
        const std::vector<double>& u, std::vector<double>& grad, void* data)
    {
        assert(data); // Need the inital displacments to compute anything.
        ObjectiveData* obj_data = static_cast<ObjectiveData*>(data);

        // Map the input displacment ot a Eigen Vector.
        const Eigen::VectorXd U
            = Eigen::Map<const Eigen::VectorXd>(u.data(), long(u.size()));
        // Compute the difference between the original and input displacments.
        const Eigen::VectorXd DIFF_U = U - obj_data->U_FLAT;
        // Only compute the gradient as needed.
        if (!grad.empty()) {
            // This should always be true according to NLopt.
            assert(size_t(grad.size()) == size_t(DIFF_U.size()));
            // ∇f(U) = ∇ (U - U0)^T * (U - U0) = 2 * (U - U0)
            Eigen::VectorXd::Map(&grad[0], long(grad.size())) = 2 * DIFF_U;

            // Save objective data for plotting later
            obj_data->objective_per_iteration.push_back(DIFF_U.squaredNorm());
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
     * @param[in] data
     */
    void volume_constraints(unsigned m, double* results, unsigned n,
        const double* u, double* grad, void* data)
    {
        assert(data); // Need data to compute volumes
        ConstraintData* constraint_data = static_cast<ConstraintData*>(data);

        // Reshape the input displacments to a |V| x 2 matrix.
        Eigen::MatrixXd U = Eigen::Map<const Eigen::VectorXd>(u, n);
        U.resize(U.rows() / 2, 2);

        // Detect the collisions for the input displacments.
        // const EdgeEdgeImpacts& ee_impacts = constraint_data->EE_IMPACTS;
        // const Eigen::VectorXi& edge_impact_map
        //     = constraint_data->EDGE_IMPACT_MAP;
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(constraint_data->E.rows());
        detect_collisions(constraint_data->V, U, constraint_data->E,
            constraint_data->CCD_DETECTION_METHOD, ee_impacts, edge_impact_map);

        // Compute the STIVs
        Eigen::VectorXd volumes;
        ccd::autodiff::compute_volumes_refresh_toi(constraint_data->V, U,
            constraint_data->E, ee_impacts, edge_impact_map,
            constraint_data->VOLUME_EPSILON, volumes);

        // Only compute the gradient as needed.
        if (grad) {
            // Gradient of the volumes (n x m)
            Eigen::MatrixXd volume_gradient;
            ccd::autodiff::compute_volumes_gradient(constraint_data->V, U,
                constraint_data->E, ee_impacts, edge_impact_map,
                constraint_data->VOLUME_EPSILON, volume_gradient);
            // This should always be true according to NLopt
            assert(size_t(m * n) == size_t(volume_gradient.size()));
            volume_gradient.resize(m * n, 1);
            Eigen::VectorXd::Map(grad, m * n) = -volume_gradient;

            // Save constraint data for plotting later
            constraint_data->sum_constraints_per_iteration.push_back(
                -volumes.sum());
        }

        // We want volumes[i] ≥ 0, but NLopt expects constraints(x) ≤ 0,
        // so negate the volume and volume gradient.
        assert(size_t(volumes.size()) == size_t(m));
        // Store the volumes in the results.
        Eigen::VectorXd::Map(results, m) = -volumes;
    }

    // Optimize the displacments using NLopt.
    bool solve_problem_with_nlopt(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const nlopt::algorithm opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt)
    {
        Eigen::MatrixXd U_flat = U;
        U_flat.resize(U.size(), 1);

        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, ccd_detection_method, ee_impacts, edge_impact_map);

        std::vector<double> objective_per_iteration;
        // Construct the data class to pass to the objective.
        ObjectiveData objectve_data(U_flat, objective_per_iteration);

        std::vector<double> sum_constraints_per_iteration;
        // Construct the data class to pass to the constraints.
        ConstraintData constraint_data(V, E, volume_epsilon,
            ccd_detection_method, ee_impacts, edge_impact_map,
            sum_constraints_per_iteration);

        nlopt::opt opt(opt_method, unsigned(U.size())); // Optimization object
        opt.set_min_objective(objective, &objectve_data);
        std::vector<double> tol(unsigned(E.rows()), 1e-8); // Tolerances

        // m volume constraints
        opt.add_inequality_mconstraint(
            volume_constraints, &constraint_data, tol);

        // Stopping criteria
        double epsilon = 1e-8;
        opt.set_xtol_rel(epsilon);
        opt.set_maxeval(int(max_iter));
        // opt.set_maxtime(1);

        // Initial guess is the value of Uopt
        std::vector<double> u(Uopt.data(), Uopt.data() + Uopt.size());
        double minf; // The minimum objective value, upon return

        // Optimize the displacments
        nlopt::result result = opt.optimize(u, minf);
        print_nlopt_termination_reason(opt, result);

        // Save JSON file of optimization objectives per iteration.
        assert(objective_per_iteration.size()
            == sum_constraints_per_iteration.size());
        OptimizationMethod m_method = opt_method == nlopt::LD_MMA ? MMA : SLSQP;
        export_intermediate(
            m_method, objective_per_iteration, sum_constraints_per_iteration);

        // Store solution in Uopt
        Uopt.col(0) = Eigen::Map<Eigen::VectorXd>(u.data(), Uopt.rows());
        Uopt.col(1) = Eigen::Map<Eigen::VectorXd>(
            &(u.data()[Uopt.rows()]), Uopt.rows());

        // Recompute the collision volumes to make sure the constraints were
        // satisfied.
        detect_collisions(
            V, Uopt, E, ccd_detection_method, ee_impacts, edge_impact_map);
        Eigen::VectorXd volumes;
        ccd::compute_volumes_fixed_toi(
            V, Uopt, E, ee_impacts, edge_impact_map, volume_epsilon, volumes);

        return result > 0 && minf >= 0
            && volumes.cwiseAbs().maxCoeff() < epsilon;
    }

    void print_nlopt_termination_reason(
        const nlopt::opt& opt, const nlopt::result result)
    {
        switch (result) {
        case nlopt::XTOL_REACHED:
            std::cout << "Optimization terminated because the relative change "
                         "was less than "
                      << opt.get_xtol_rel() << "." << std::endl;
            break;
        case nlopt::MAXEVAL_REACHED:
            std::cout << "Optimization terminated because the number of "
                         "evaluations exceeded "
                      << opt.get_maxeval() << "." << std::endl;
            break;
        case nlopt::MAXTIME_REACHED:
            std::cout << "Optimization terminated because the runtime exceeded "
                      << opt.get_maxtime() << " seconds." << std::endl;
            break;
        default:
            std::cout << "Optimization terminated because of code " << result
                      << " (see "
                         "https://nlopt.readthedocs.io/en/latest/"
                         "NLopt_Reference/#return-values)."
                      << std::endl;
            break;
        }
    }

} // namespace opt
} // namespace ccd
