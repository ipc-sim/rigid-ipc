#include "sympletic_euler_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    // Take a single time step in 2D.
    void SympleticEulerTimeStepper::step2D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector2d& gravity,
        const double& time_step) const
    {
        physics::Pose<double> zero = physics::Pose<double>::Zero(bodies.dim());
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Store the previous configurations and velocities
            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // Integrate linear terms
            rb.velocity.position = rb.velocity_prev.position
                + time_step * (gravity + rb.force.position / rb.mass);
            // Fix linear terms
            rb.velocity.position =
                (rb.is_dof_fixed.head(rb.velocity.pos_ndof()))
                    .select(zero.position, rb.velocity.position);
            // Update position using update velocity
            rb.pose.position =
                rb.pose_prev.position + time_step * rb.velocity.position;

            // Integrate angular terms
            rb.velocity.rotation = rb.velocity_prev.rotation
                + time_step
                    * rb.moment_of_inertia.lu().solve(rb.force.rotation);
            // Fix angular terms
            rb.velocity.rotation =
                (rb.is_dof_fixed.tail(rb.velocity.rot_ndof()))
                    .select(zero.rotation, rb.velocity.rotation);
            // Update position using update velocity
            rb.pose.rotation =
                rb.pose_prev.rotation + time_step * rb.velocity.rotation;
        });
    }

} // namespace time_stepper
} // namespace ccd
