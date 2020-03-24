#include "verlet_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    // Take a single time step in 2D.
    void VerletTimeStepper::step2D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector2d& gravity,
        const double& time_step) const
    {
        physics::Pose<double> zero = physics::Pose<double>::Zero(bodies.dim());
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Zero out velocity of fixed dof
            rb.velocity.position =
                (rb.is_dof_fixed.head(rb.velocity.pos_ndof()))
                    .select(zero.position, rb.velocity.position);
            rb.velocity.rotation =
                (rb.is_dof_fixed.head(rb.velocity.rot_ndof()))
                    .select(zero.rotation, rb.velocity.rotation);

            // Store the previous configurations and velocities
            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // F( q0 )
            // fsys.computeForce(q0, v0, start_time, F);
            physics::Pose<double> F = rb.force;
            F.position += rb.mass * gravity;
            F.position = (rb.is_dof_fixed.head(F.pos_ndof()))
                             .select(zero.position, F.position);
            F.rotation = (rb.is_dof_fixed.head(F.rot_ndof()))
                             .select(zero.rotation, F.rotation);

            // v_{1/2} = v0 + (1/2) h a0
            rb.velocity += 0.5 * time_step * F / rb.mass;

            // q1 = q0 + h v0 + (1/2) h^2 a0
            rb.pose += time_step * rb.velocity;

            // F( q1 )
            // fsys.computeForce(q1, v0, start_time, F);
            // zeroForcesOnKinematicBodies(fsys, F);

            // v1 = v0 + (1/2) h ( a0 + a1 )
            rb.velocity += 0.5 * time_step * F / rb.mass;
        });
    }

} // namespace time_stepper
} // namespace ccd
