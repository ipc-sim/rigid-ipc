#include "sympletic_euler_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    // Take a single time step in 2D.
    void SympleticEulerTimeStepper::step2D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector2d& gravity,
        const double& time_step) const
    {
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Store the previous configurations and velocities
            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // Update linear velocities
            rb.velocity.position +=
                time_step * (gravity + rb.force.position / rb.mass);
            // Update angular velocity
            rb.velocity.rotation += time_step
                * rb.moment_of_inertia.cwiseInverse().asDiagonal()
                * rb.force.rotation;
            // Zero out velocity of fixed dof
            rb.velocity.zero_dof(rb.is_dof_fixed);

            // Update position using update velocity
            rb.pose.position += time_step * rb.velocity.position;
            // Update position using update velocity
            rb.pose.rotation += time_step * rb.velocity.rotation;
        });
    }

} // namespace time_stepper
} // namespace ccd
