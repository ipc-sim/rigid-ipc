#include "sympletic_euler_time_stepper.hpp"

namespace ipc::rigid {

// Take a single time step in 2D.
void SympleticEulerTimeStepper::step2D(
    RigidBody& body,
    const Eigen::Vector2d& gravity,
    const double& time_step) const
{
    assert(body.dim() == 2);
    // Store the previous configurations and velocities
    body.pose_prev = body.pose;
    body.velocity_prev = body.velocity;

    // Update linear velocities
    body.velocity.position +=
        time_step * (gravity + body.force.position / body.mass);
    // Update angular velocity
    body.velocity.rotation += time_step
        * body.moment_of_inertia.cwiseInverse().asDiagonal()
        * body.force.rotation;
    // Zero out velocity of fixed dof
    body.velocity.zero_dof(body.is_dof_fixed);

    // Update position using update velocity
    body.pose.position += time_step * body.velocity.position;
    // Update position using update velocity
    body.pose.rotation += time_step * body.velocity.rotation;
}

} // namespace ipc::rigid
