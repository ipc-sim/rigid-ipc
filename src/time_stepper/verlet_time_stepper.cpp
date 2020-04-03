#include "verlet_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    // Take a single time step in 2D.
    void VerletTimeStepper::step2D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector2d& gravity,
        const double& time_step) const
    {
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Zero out velocity of fixed dof
            rb.velocity.zero_dof(rb.is_dof_fixed);

            // Store the previous configurations and velocities
            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // F( q0 )
            // fsys.computeForce(q0, v0, start_time, F);
            physics::Pose<double> F = rb.force;
            F.position += rb.mass * gravity;
            // Zero out forces of fixed dof
            F.zero_dof(rb.is_dof_fixed);

            // v_{1/2} = v0 + (1/2) h a0
            rb.velocity.position += 0.5 * time_step * F.position / rb.mass;
            rb.velocity.rotation += 0.5 * time_step * F.rotation / rb.mass;

            // q1 = q0 + h v0 + (1/2) h^2 a0
            rb.pose.position += time_step * rb.velocity.position;
            rb.pose.rotation += time_step * rb.velocity.rotation;

            // F( q1 )
            // fsys.computeForce(q1, v0, start_time, F);
            // zeroForcesOnKinematicBodies(fsys, F);

            // v1 = v0 + (1/2) h ( a0 + a1 )
            rb.velocity.position += 0.5 * time_step * F.position / rb.mass;
            rb.velocity.rotation += 0.5 * time_step * F.rotation / rb.mass;
        });
    }

} // namespace time_stepper
} // namespace ccd
