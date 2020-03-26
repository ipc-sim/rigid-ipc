#include "exponential_euler_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    static void project_orientation(Eigen::Matrix3d& R)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> svd;
        svd.compute(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
        R = svd.matrixU() * svd.matrixV().transpose();
        assert(R.isUnitary(1e-9));
        assert(fabs(R.determinant() - 1.0) < 1.0e-6);
    }

    void ExponentialEulerTimeStepper::step3D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector3d& gravity,
        const double& time_step) const
    {
        assert(bodies.dim() == 3);

        // Store the previous configurations and velocities
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Zero out velocity of fixed dof
            // Fixed dof is specified in body frame so provide a transform
            rb.velocity.zero_dof(rb.is_dof_fixed, rb.R0);

            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // Update the center of mass position
            rb.pose.position += time_step * rb.velocity.position;

            // Update the orientaiton
            // R₁ = R₀ + h * R₀ω̂
            Eigen::Matrix3d R1 = rb.pose_prev.construct_rotation_matrix();
            R1 += time_step * R1 * Eigen::Hat(rb.velocity.rotation);
            project_orientation(R1);
            rb.pose.rotation = R1.eulerAngles(2, 1, 0).reverse();

            // Compute the acceleration at ( q0, v0 )
            physics::Pose<double> acceleration(
                rb.force.position / rb.mass,
                rb.moment_of_inertia.cwiseInverse().asDiagonal()
                    * rb.force.rotation);
            acceleration.position += gravity;
            // Zero out acceleration of fixed dof
            // Fixed dof is specified in body frame so provide a transform
            acceleration.zero_dof(rb.is_dof_fixed, rb.R0);

            // Update the velocity
            rb.velocity += time_step * acceleration;
        });
    }

} // namespace time_stepper
} // namespace ccd
