#pragma once

#include <autodiff/autodiff_types.hpp>
#include <multiprecision.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <physics/rigid_body_problem.hpp>
#include <solvers/barrier_solver.hpp>

namespace ccd {

namespace opt {

    struct RigidBodyEdgeVertexCandidate {
        long vertex_body_id;
        long edge_body_id;
        long vertex_local_id;
        long edge_vertex0_local_id;
        long edge_vertex1_local_id;

        inline std::array<long, 2> get_body_ids()
        {
            return { { vertex_body_id, edge_body_id } };
        }
    };
    struct RigidBodyEdgeEdgeCandidate {
        long edge0_body_id;
        long edge1_body_id;
        long edge0_vertex0_local_id;
        long edge0_vertex1_local_id;
        long edge1_vertex0_local_id;
        long edge1_vertex1_local_id;

        inline std::array<long, 2> get_body_ids()
        {
            return { { edge0_body_id, edge1_body_id } };
        }
    };
    struct RigidBodyFaceVertexCandidate {
        long vertex_body_id;
        long face_body_id;
        long vertex_local_id;
        long face_vertex0_local_id;
        long face_vertex1_local_id;
        long face_vertex2_local_id;

        inline std::array<long, 2> get_body_ids()
        {
            return { { vertex_body_id, face_body_id } };
        }
    };

    /// This class is both a simulation and optimization problem.
    class DistanceBarrierRBProblem : public physics::RigidBodyProblem,
                                     public virtual BarrierProblem {
    public:
        DistanceBarrierRBProblem();
        virtual ~DistanceBarrierRBProblem() = default;

        void settings(const nlohmann::json& params) override;
        nlohmann::json state() const override;

        static std::string problem_name()
        {
            return "distance_barrier_rb_problem";
        }

        virtual std::string name() const override
        {
            return DistanceBarrierRBProblem::problem_name();
        }

        ////////////////////////////////////////////////////////////
        /// Rigid Body Problem

        bool simulation_step(const double time_step) override;
        bool take_step(
            const Eigen::VectorXd& sigma, const double time_step) override;

        ////////////////////////////////////////////////////////////
        /// Optimization Problem

        /// @returns the number of variables
        int num_vars() const override { return num_vars_; }

        /// @returns A vector of booleans indicating if a DoF is fixed.
        const Eigen::VectorXb& is_dof_fixed() override
        {
            return m_assembler.is_rb_dof_fixed;
        }

        /// Determine if there is a collision between two configurations
        bool has_collisions(
            const Eigen::VectorXd& sigma_i,
            const Eigen::VectorXd& sigma_j) const override;

        ////////////////////////////////////////////////////////////
        /// Barrier Problem

        /// Compute E(x) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        double compute_energy_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad = true,
            bool compute_hess = true) override;

        /// Compute ∑_{k ∈ C} b(d(x_k)) in f(x) = E(x) + κ ∑_{k ∈ C} b(d(x_k))
        double compute_barrier_term(
            const Eigen::VectorXd& x,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            int& num_constraints,
            bool compute_grad = true,
            bool compute_hess = true) override;

        // Include thes lines to avoid issues with overriding inherited
        // functions with the same name.
        // (http://www.cplusplus.com/forum/beginner/24978/)
        using BarrierProblem::compute_barrier_term;
        using BarrierProblem::compute_energy_term;

        double
        compute_min_distance(const Eigen::VectorXd& sigma) const override;

        double get_barrier_homotopy() const override
        {
            return constraint_.get_barrier_epsilon();
        }
        void set_barrier_homotopy(const double eps) override
        {
            constraint_.set_barrier_epsilon(eps);
        }

        double get_barrier_stiffness() const override
        {
            return barrier_stiffness;
        }
        void set_barrier_stiffness(const double kappa) override
        {
            barrier_stiffness = kappa;
        }

        opt::CollisionConstraint& constraint() override { return constraint_; }
        const opt::CollisionConstraint& constraint() const override
        {
            return constraint_;
        }
        opt::OptimizationSolver& solver() override { return opt_solver_; }

    protected:
        void extract_local_system(
            const EdgeVertexCandidate& c, RigidBodyEdgeVertexCandidate& rbc);
        void extract_local_system(
            const EdgeEdgeCandidate& c, RigidBodyEdgeEdgeCandidate& rbc);
        void extract_local_system(
            const FaceVertexCandidate& c, RigidBodyFaceVertexCandidate& rbc);

        template <typename T, typename RigidBodyCandidate>
        T distance_barrier(
            const Eigen::VectorXd& sigma, const RigidBodyCandidate& rbc);

        template <typename T>
        T distance(
            const Eigen::VectorXd& sigma,
            const RigidBodyEdgeVertexCandidate& rbc);
        template <typename T>
        T distance(
            const Eigen::VectorXd& sigma,
            const RigidBodyEdgeEdgeCandidate& rbc);
        template <typename T>
        T distance(
            const Eigen::VectorXd& sigma,
            const RigidBodyFaceVertexCandidate& rbc);

        template <typename Candidate, typename RigidBodyCandidate>
        void add_constraint_barrier(
            const Eigen::VectorXd& sigma,
            const Candidate& candidate,
            double& Bx,
            Eigen::VectorXd& grad,
            std::vector<Eigen::Triplet<double>>& hess_triplets,
            bool compute_grad,
            bool compute_hess);

        /// Computes the barrier term value, gradient, and hessian
        double compute_barrier_term(
            const Eigen::VectorXd& sigma,
            const Candidates& distance_candidates,
            Eigen::VectorXd& grad,
            Eigen::SparseMatrix<double>& hess,
            bool compute_grad,
            bool compute_hess);

#ifdef WITH_DERIVATIVE_CHECK
        // The following functions are used exclusivly to check that the
        // gradient and hessian match a finite difference version.

        Eigen::VectorXd compute_full_grad_barrier(
            const Eigen::VectorXd& sigma, const Candidates& candidates);

        template <typename Candidate, typename RigidBodyCandidate>
        void check_distance_finite_diff(
            const Eigen::VectorXd& sigma, const Candidate& candidate);

        void check_grad_barrier(
            const Eigen::VectorXd& sigma,
            const Candidates& candidates,
            const Eigen::VectorXd& grad);
#endif

        double min_distance;
        opt::DistanceBarrierConstraint constraint_;
        opt::BarrierSolver opt_solver_;

        // double barrier_homotopy; ///< \f$\hat{d}\f$
        double barrier_stiffness; ///< \f$\kappa\f$
    };

} // namespace opt
} // namespace ccd
