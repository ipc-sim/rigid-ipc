#pragma once

#include <time_stepper/time_stepper.hpp>

namespace ccd {
namespace time_stepper {

    class DMVTimeStepper : public TimeStepper {
    public:
        virtual ~DMVTimeStepper() = default;

        static std::string default_name() { return "dmv"; }
        virtual std::string name() const override
        {
            return DMVTimeStepper::default_name();
        }

    protected:
        /**
         * @brief Take a single time step.
         *
         * @param bodies     Rigid bodies
         * @param gravity    Acceleration due to gravity
         * @param time_step  Timestep
         */
        virtual void step3D(
            physics::RigidBodyAssembler& bodies,
            const Eigen::Vector3d& gravity,
            const double& time_step) const override;
    };

} // namespace time_stepper
} // namespace ccd
