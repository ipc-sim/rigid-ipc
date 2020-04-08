#include "dmv_time_stepper.hpp"

#include <logger.hpp>

namespace ccd {
namespace time_stepper {

    void solve_DMV(
        const Eigen::Vector3d& angular_momentum,
        const double& time_step,
        const Eigen::Vector3d& inertia,
        Eigen::Quaterniond& q)
    {
        const double eps = fabs(time_step * 1e-15);
        const double ha = time_step / 2.0;
        const double fac1 = (inertia.y() - inertia.z()) / inertia.x();
        const double fac2 = (inertia.z() - inertia.x()) / inertia.y();
        const double fac3 = (inertia.x() - inertia.y()) / inertia.z();

        const double am1i = angular_momentum.x() * ha / inertia.x();
        const double am2i = angular_momentum.y() * ha / inertia.y();
        const double am3i = angular_momentum.z() * ha / inertia.z();
        double cm1 = am1i + fac1 * am2i * am3i;
        double cm2 = am2i + fac2 * cm1 * am3i;
        double cm3 = am3i + fac3 * cm1 * cm2;
        double err = std::numeric_limits<double>::infinity();
        for (unsigned dmv_itr = 0; dmv_itr < 50; ++dmv_itr) {
            const double cm1b = cm1;
            const double cm2b = cm2;
            const double cm3b = cm3;
            const double calpha = cm1 * cm1 + 1.0 + cm2 * cm2 + cm3 * cm3;
            cm1 = calpha * am1i + fac1 * cm2 * cm3;
            cm2 = calpha * am2i + fac2 * cm1 * cm3;
            cm3 = calpha * am3i + fac3 * cm1 * cm2;
            err = fabs(cm1b - cm1) + fabs(cm2b - cm2) + fabs(cm3b - cm3);
            if (err <= eps) {
                break;
            }
        }
        if (err > eps || !isfinite(err)) {
            spdlog::warn("DMV failed to terminate");
        }
        assert(isfinite(err));

        const double q0 = q.w();
        const double q1 = q.x();
        const double q2 = q.y();
        const double q3 = q.z();
        q.w() = q0 - cm1 * q1 - cm2 * q2 - cm3 * q3;
        q.x() = q1 + cm1 * q0 + cm3 * q2 - cm2 * q3;
        q.y() = q2 + cm2 * q0 + cm1 * q3 - cm3 * q1;
        q.z() = q3 + cm3 * q0 + cm2 * q1 - cm1 * q2;

        q.normalize();
    }

#ifndef NDEBUG
    bool is_unit_length(const Eigen::Quaterniond& q)
    {
        const double squared_norm =
            q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z();
        return fabs(squared_norm - 1.0) <= 1.0e-9;
    }
#endif

    void
    DMV(const Eigen::Matrix3d& R0,
        const Eigen::Vector3d& angular_momentum,
        double time_step,
        const Eigen::Vector3d& moment_of_inertia,
        Eigen::Matrix3d& R1)
    {
        // Get the orientation of the body
        assert(fabs(R0.determinant() - 1.0) <= 1.0e-9);
        assert(R0.isUnitary(1e-9));

        // Convert the input orientation to a quaternion
        Eigen::Quaterniond Q_new(R0);
        assert(is_unit_length(Q_new));

        // Solve for the next orientation
        solve_DMV(angular_momentum, time_step, moment_of_inertia, Q_new);
        assert(is_unit_length(Q_new));

        // Convert the output orientation to a matrix
        R1 = Q_new.toRotationMatrix();
        // Ensure we have still have an orthonormal rotation matrix
        assert(fabs(R1.determinant() - 1.0) <= 1.0e-9);
        assert(R1.isUnitary(1e-9));
    }

    /**
     * @brief Take a single time step.
     *
     * @param bodies     Rigid bodies
     * @param gravity    Acceleration due to gravity
     * @param time_step  Timestep
     */
    void DMVTimeStepper::step3D(
        physics::RigidBodyAssembler& bodies,
        const Eigen::Vector3d& gravity,
        const double& time_step) const
    {
        assert(bodies.dim() == 3);
        physics::Pose<double> zero = physics::Pose<double>::Zero(bodies.dim());

        // Store the previous configurations and velocities
        tbb::parallel_for_each(bodies.m_rbs, [&](physics::RigidBody& rb) {
            // Zero out velocity of fixed dof
            // Fixed dof is specified in body frame so provide a transform
            rb.velocity.zero_dof(rb.is_dof_fixed, rb.R0);

            rb.pose_prev = rb.pose;
            rb.velocity_prev = rb.velocity;

            // p1 = m * v0
            physics::Pose<double> momentum_t1(
                rb.mass * rb.velocity.position,
                rb.moment_of_inertia.asDiagonal() * rb.velocity.rotation);

            // TODO: Compute forces at start of time_step
            // Hamiltonian so there shouldn't be velocity dependent forces
            // fsys.computeForce(q0, v0, next_time, F);
            physics::Pose<double> F = rb.force;
            F.position += rb.mass * gravity; // mass will be divided out
            // Zero out force of fixed dof
            // Fixed dof is specified in body frame so provide a transform
            F.zero_dof(rb.is_dof_fixed, rb.R0);

            // First momentum update
            // p1 += 0.5 * h * F_q0;
            momentum_t1.position += 0.5 * time_step * F.position;
            momentum_t1.rotation += 0.5 * time_step * F.rotation;

            // Linear position update
            rb.pose.position += time_step * momentum_t1.position / rb.mass;

            // Angular position update
            Eigen::Matrix3d R1;
            DMV(rb.pose.construct_rotation_matrix(), momentum_t1.rotation,
                time_step, rb.moment_of_inertia, R1);
            Eigen::AngleAxisd r1 = Eigen::AngleAxisd(R1);
            rb.pose.rotation = r1.angle() * r1.axis();

            // TODO: Compute forces at end of time_step
            // Hamiltonian so there shouldn't be velocity dependent forces
            // fsys.computeForce(q1, v0, next_time, F);

            // Second momentum update
            // p1 += 0.5 * h * F_q1;
            momentum_t1.position += 0.5 * time_step * F.position;
            momentum_t1.rotation += 0.5 * time_step * F.rotation;

            // Convert momentum to velocity
            // Linear component
            rb.velocity.position = momentum_t1.position / rb.mass;
            // Rotational component
            Eigen::Matrix3d R = rb.pose.construct_rotation_matrix();
            assert(fabs(R.determinant() - 1.0) < 1.0e-9);
            assert(R.isUnitary(1e-9));
            Eigen::Vector3d Iinv0 = rb.moment_of_inertia.cwiseInverse();
            assert((Iinv0.array() > 0.0).all());
            rb.velocity.rotation = Iinv0.asDiagonal() * momentum_t1.rotation;
        });
    }

} // namespace time_stepper
} // namespace ccd
