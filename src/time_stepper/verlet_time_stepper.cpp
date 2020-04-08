#include "verlet_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    // Take a single time step in 2D.
    void VerletTimeStepper::step2D(
        physics::RigidBody& body,
        const Eigen::Vector2d& gravity,
        const double& time_step) const
    {
        assert(body.dim() == 2);

        // Zero out velocity of fixed dof
        body.velocity.zero_dof(body.is_dof_fixed);

        // Store the previous configurations and velocities
        body.pose_prev = body.pose;
        body.velocity_prev = body.velocity;

        // F( q0 )
        // fsys.computeForce(q0, v0, start_time, F);
        physics::Pose<double> F = body.force;
        F.position += body.mass * gravity;
        // Zero out forces of fixed dof
        F.zero_dof(body.is_dof_fixed);

        // v_{1/2} = v0 + (1/2) h a0
        body.velocity.position += 0.5 * time_step * F.position / body.mass;
        body.velocity.rotation += 0.5 * time_step * F.rotation / body.mass;

        // q1 = q0 + h v0 + (1/2) h^2 a0
        body.pose.position += time_step * body.velocity.position;
        body.pose.rotation += time_step * body.velocity.rotation;

        // F( q1 )
        // fsys.computeForce(q1, v0, start_time, F);
        // zeroForcesOnKinematicBodies(fsys, F);

        // v1 = v0 + (1/2) h ( a0 + a1 )
        body.velocity.position += 0.5 * time_step * F.position / body.mass;
        body.velocity.rotation += 0.5 * time_step * F.rotation / body.mass;
    }

} // namespace time_stepper
} // namespace ccd
