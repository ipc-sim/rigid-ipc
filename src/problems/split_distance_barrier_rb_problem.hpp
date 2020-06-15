#pragma once

#include <problems/distance_barrier_rb_problem.hpp>

namespace ccd {
namespace opt {

    /// @brief A diststance barrier rigid body problem using a split time
    /// stepping method.
    class SplitDistanceBarrierRBProblem : public DistanceBarrierRBProblem {
    public:
        SplitDistanceBarrierRBProblem();
        virtual ~SplitDistanceBarrierRBProblem() = default;

        /// The name of the class as a string
        static std::string problem_name()
        {
            return "split_distance_barrier_rb_problem";
        }

        /// The name of the class as a string
        virtual std::string name() const override
        {
            return SplitDistanceBarrierRBProblem::problem_name();
        }

        ////////////////////////////////////////////////////////////
        // Rigid Body Problem

        void simulation_step(
            bool& had_collisions,
            bool& has_intersections,
            bool solve_collisions = true) override;

        ////////////////////////////////////////////////////////////
        // Barrier Problem

        /// @brief Compute \f$E(x)\f$ in
        /// \f$f(x) = E(x) + \kappa \sum_{k \in C} b(d(x_k))\f$
        double compute_energy_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad = true,
            bool compute_hess = true) override;

        // Include thes lines to avoid issues with overriding inherited
        // functions with the same name.
        // (http://www.cplusplus.com/forum/beginner/24978/)
        using BarrierProblem::compute_barrier_term;
        using BarrierProblem::compute_energy_term;
    };

} // namespace opt
} // namespace ccd
