#include "exponential_euler_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

    static void project_orientation(Eigen::Matrix3d& R)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> svd;
        svd.compute(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
        R = svd.matrixU() * svd.matrixV().transpose();
        // clang-format off
        assert((R * R.transpose() - Eigen::Matrix3d::Identity())
                .lpNorm<Eigen::Infinity>() < 1.0e-9);
        // clang-format on
        assert(fabs(R.determinant() - 1.0) < 1.0e-6);
    }

    void ExponentialEulerTimeStepper::step3D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector3d& gravity,
        const double& time_step) const
    {
        assert(bodies.dim() == 3);
        physics::Pose<double> zero = physics::Pose<double>::Zero(bodies.dim());

        // Store the previous configurations and velocities
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // Update the center of mass position
            rb.pose.position += time_step * rb.velocity.position;

            // Update the orientaiton
            // R₁ = R₀ + h * ω̂R₀
            Eigen::Matrix3d R1 = rb.pose_prev.construct_rotation_matrix();
            R1 += time_step * Eigen::Hat(rb.velocity.rotation) * R1;
            project_orientation(R1);
            rb.pose.rotation = R1.eulerAngles(2, 1, 0);

            // Compute the acceleration at ( q0, v0 )
            physics::Pose<double> acceleration(
                rb.force.position / rb.mass,
                rb.moment_of_inertia.cwiseInverse().asDiagonal()
                    * rb.force.rotation);
            acceleration.position += gravity;

            // Update the velocity
            rb.velocity += acceleration * time_step;
        });
    }

} // namespace time_stepper
} // namespace ccd
