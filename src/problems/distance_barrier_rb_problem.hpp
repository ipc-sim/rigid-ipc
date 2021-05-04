#pragma once

#include <tbb/concurrent_vector.h>

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/friction_constraint.hpp>

#include <autodiff/autodiff_types.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/optimization_problem.hpp>
#include <physics/rigid_body_problem.hpp>
#include <problems/rigid_body_collision_constraint.hpp>
#include <solvers/homotopy_solver.hpp>
#include <utils/multiprecision.hpp>

namespace ipc::rigid {

/// @brief Possible methods for detecting all edge vertex collisions.
enum BodyEnergyIntegrationMethod {
    IMPLICIT_EULER,
    IMPLICIT_NEWMARK,
    STABILIZED_NEWMARK
};

const static BodyEnergyIntegrationMethod
    DEFAULT_BODY_ENERGY_INTEGRATION_METHOD = IMPLICIT_EULER;

NLOHMANN_JSON_SERIALIZE_ENUM(
    BodyEnergyIntegrationMethod,
    { { IMPLICIT_EULER, "implicit_euler" },
      { IMPLICIT_NEWMARK, "implicit_newmark" },
      { STABILIZED_NEWMARK, "stabilized_newmark" },
      { DEFAULT_BODY_ENERGY_INTEGRATION_METHOD, "default" } });

/// This class is both a simulation and optimization problem.
class DistanceBarrierRBProblem : public RigidBodyProblem,
                                 public virtual BarrierProblem {
public:
    DistanceBarrierRBProblem();
    virtual ~DistanceBarrierRBProblem() = default;

    bool settings(const nlohmann::json& params) override;
    nlohmann::json settings() const override;

    nlohmann::json state() const override;

    static std::string problem_name() { return "distance_barrier_rb_problem"; }

    virtual std::string name() const override
    {
        return DistanceBarrierRBProblem::problem_name();
    }

    ////////////////////////////////////////////////////////////
    // Rigid Body Problem

    void simulation_step(
        bool& had_collisions,
        bool& has_intersections,
        bool solve_collisions = true) override;

    bool take_step(const Eigen::VectorXd& x) override;

    /// Use the solver to solve this problem.
    OptimizationResults solve_constraints() override;

    ////////////////////////////////////////////////////////////
    // Optimization Problem

    /// @returns the number of variables
    int num_vars() const override { return num_vars_; }

    /// @returns A vector of booleans indicating if a DoF is fixed.
    const VectorXb& is_dof_fixed() const override
    {
        return m_assembler.is_rb_dof_fixed;
    }

    Eigen::VectorXi free_dof() const override;

    /// Determine if there is a collision between two configurations
    bool has_collisions(
        const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j) override;

    /// Compute the earliest time of impact between two configurations
    double compute_earliest_toi(
        const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j) override;

    bool is_ccd_aligned_with_newton_update() override
    {
        return m_constraint.trajectory_type != TrajectoryType::LINEAR;
    }

    /// Get the world coordinates of the vertices
    Eigen::MatrixXd world_vertices(const Eigen::VectorXd& x) const override
    {
        return m_assembler.world_vertices(this->dofs_to_poses(x));
    }

    /// Get the length of the diagonal of the worlds bounding box
    double world_bbox_diagonal() const override
    {
        // TODO: Compute this value dynamicly if necessary
        return init_bbox_diagonal;
    }

    /// Get the mass matrix
    DiagonalMatrixXd mass_matrix() const override
    {
        return m_assembler.m_rb_mass_matrix;
    }

    /// Get the average mass (average of mass matrix diagonal)
    double average_mass() const override { return m_assembler.average_mass; }

    /// Get the time-step
    double timestep() const override { return RigidBodyProblem::timestep(); }

    ////////////////////////////////////////////////////////////
    // Barrier Problem

    /// Compute the objective function f(x)
    double compute_objective(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true) override;

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

    double compute_friction_term(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad = true,
        bool compute_hess = true);

    virtual double compute_friction_term(const Eigen::VectorXd& x) final
    {
        Eigen::VectorXd grad;
        Eigen::SparseMatrix<double> hess;
        return compute_friction_term(
            x, grad, hess, /*compute_grad=*/false, /*compute_hess=*/false);
    }

    virtual double
    compute_friction_term(const Eigen::VectorXd& x, Eigen::VectorXd& grad) final
    {
        Eigen::SparseMatrix<double> hess;
        return compute_friction_term(
            x, grad, hess, /*compute_grad=*/true, /*compute_hess=*/false);
    }

    // Include thes lines to avoid issues with overriding inherited
    // functions with the same name.
    // (http://www.cplusplus.com/forum/beginner/24978/)
    using BarrierProblem::compute_barrier_term;
    using BarrierProblem::compute_energy_term;

    /// Compute the minimum distance among geometry
    double compute_min_distance() const override;
    double compute_min_distance(const Eigen::VectorXd& x) const override;

    /// Compute the value of the barrier at a distance x
    double barrier_hessian(double x) const override
    {
        return m_constraint.distance_barrier_hessian(x);
    }

    double barrier_activation_distance() const override
    {
        return m_constraint.barrier_activation_distance();
    }
    void barrier_activation_distance(double dhat) override
    {
        m_constraint.barrier_activation_distance(dhat);
    }

    double barrier_stiffness() const override { return m_barrier_stiffness; }
    void barrier_stiffness(double kappa) override
    {
        m_barrier_stiffness = kappa;
    }

    CollisionConstraint& constraint() override { return m_constraint; }
    const CollisionConstraint& constraint() const override
    {
        return m_constraint;
    }
    OptimizationSolver& solver() override { return *m_opt_solver; }

    int num_contacts() const override { return m_num_contacts; };

    ////////////////////////////////////////////////////////////
    // Augmented Lagrangian for equality constraints

    /// Initialize the augmented Lagrangian variables.
    void init_augmented_lagrangian();
    /// Update the target poses of kinematic bodies
    void step_kinematic_bodies();

    /// Update the augmented Lagrangian for kinematic bodies.
    void update_augmented_lagrangian(const Eigen::VectorXd& x) override;

    /// Compute the convergence criteria η of the augment Lagrangian.
    double
    compute_linear_augment_lagrangian_progress(const Eigen::VectorXd& x) const;
    double
    compute_angular_augment_lagrangian_progress(const Eigen::VectorXd& x) const;

    /// Determine if the value x satisfies the equality constraints.
    bool
    are_equality_constraints_satisfied(const Eigen::VectorXd& x) const override;

    /// Compute the augmented Lagrangian potential used to enforce equality
    /// constraints.
    double compute_augmented_lagrangian(
        const Eigen::VectorXd& x,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess);

protected:
    /// Update problem using current status of bodies.
    virtual void update_constraints() override;

    /// Update problem using current status of bodies.
    void update_friction_constraints(
        const Constraints& collision_constraints, const PosesD& poses);

    template <typename T>
    T compute_body_energy(
        const RigidBody& body,
        const Pose<T>& pose,
        const VectorMax6d& grad_barrier_t0);

    template <typename RigidBodyConstraint, typename FrictionConstraint>
    double compute_friction_potential(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXd& jac_V,
        const Eigen::MatrixXd& hess_V,
        const FrictionConstraint& constraint,
        Eigen::VectorXd& grad,
        std::vector<Eigen::Triplet<double>>& hess_triplets,
        bool compute_grad,
        bool compute_hess);

    /// Computes the barrier term value, gradient, and hessian from
    /// distance constraints.
    double compute_barrier_term(
        const Eigen::VectorXd& x,
        const Constraints& distance_constraints,
        Eigen::VectorXd& grad,
        Eigen::SparseMatrix<double>& hess,
        bool compute_grad,
        bool compute_hess);

    virtual double compute_barrier_term(
        const Eigen::VectorXd& x, const Constraints& distance_constraints) final
    {
        Eigen::VectorXd grad;
        Eigen::SparseMatrix<double> hess;
        return compute_barrier_term(
            x, distance_constraints, grad, hess, /*compute_grad=*/false,
            /*compute_hess=*/false);
    }

    virtual double compute_barrier_term(
        const Eigen::VectorXd& x,
        const Constraints& distance_constraints,
        Eigen::VectorXd& grad) final
    {
        Eigen::SparseMatrix<double> hess;
        return compute_barrier_term(
            x, distance_constraints, grad, hess, /*compute_grad=*/true,
            /*compute_hess=*/false);
    }

#ifdef RIGID_IPC_WITH_DERIVATIVE_CHECK
    // The following functions are used exclusivly to check that the
    // gradient and hessian match a finite difference version.

    template <typename RigidBodyConstraint, typename Constraint>
    void check_distance_finite_gradient(
        const Eigen::VectorXd& x, const Constraint& constraint);
    template <typename RigidBodyConstraint, typename Constraint>
    void check_distance_finite_hessian(
        const Eigen::VectorXd& x, const Constraint& constraint);

    void check_barrier_gradient(
        const Eigen::VectorXd& x,
        const Constraints& constraints,
        const Eigen::VectorXd& grad);
    void check_barrier_hessian(
        const Eigen::VectorXd& x,
        const Constraints& constraints,
        const Eigen::SparseMatrix<double>& hess);

    void check_friction_gradient(
        const Eigen::VectorXd& x, const Eigen::VectorXd& grad);
    void check_friction_hessian(
        const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& hess);

    void check_augmented_lagrangian_gradient(
        const Eigen::VectorXd& x, const Eigen::VectorXd& grad);
    void check_augmented_lagrangian_hessian(
        const Eigen::VectorXd& x, const Eigen::SparseMatrix<double>& hess);

    bool is_checking_derivative = false;
#endif

    /// @brief Constraint helper for active set and collision detection.
    DistanceBarrierConstraint m_constraint;

    /// @brief Solver for solving this optimization problem.
    std::shared_ptr<OptimizationSolver> m_opt_solver;

    /// @brief Multiplier of barrier term in objective, \f$\kappa\f$.
    double m_barrier_stiffness;

    /// @brief Current minimum distance between bodies.
    /// Negative values indicate a minimum distance greater than the
    /// activation distance.
    double min_distance;

    /// @brief Did the step have collisions?
    bool m_had_collisions;
    /// @brief The number of collision during the timestep.
    int m_num_contacts;

    /// @brief Gradient of barrier potential at the start of the time-step.
    Eigen::VectorXd grad_barrier_t0;

    // Friction
    double static_friction_speed_bound;
    int friction_iterations;
    FrictionConstraints friction_constraints;

    // Augmented Lagrangian
    double linear_augmented_lagrangian_penalty;
    double angular_augmented_lagrangian_penalty;
    Eigen::VectorXd linear_augmented_lagrangian_multiplier;
    Eigen::MatrixXd angular_augmented_lagrangian_multiplier;
    Eigen::VectorXd x_pred; ///< Predicted DoF using unconstrained timestep
    VectorXb is_dof_satisfied;

private:
    /// Method for integrating the body energy.
    BodyEnergyIntegrationMethod body_energy_integration_method;
};

} // namespace ipc::rigid
