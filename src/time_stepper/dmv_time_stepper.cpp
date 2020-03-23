#include "dmv_time_stepper.hpp"

namespace ccd {
namespace time_stepper {

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
        throw NotImplementedError("DMV time-stepper not implemented!");
    }

} // namespace time_stepper
} // namespace ccd
