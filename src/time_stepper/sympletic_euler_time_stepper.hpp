#pragma once

#include <time_stepper/time_stepper.hpp>

namespace ccd {
namespace time_stepper {

    class SympleticEulerTimeStepper : public TimeStepper {
    public:
        virtual ~SympleticEulerTimeStepper() = default;

        static std::string default_name() { return "sympletic_euler"; }
        virtual std::string name() const override
        {
            return SympleticEulerTimeStepper::default_name();
        }

    protected:
        /**
         * @brief Take a single time step.
         *
         * @param bodies     Rigid bodies
         * @param gravity    Acceleration due to gravity
         * @param time_step  Timestep
         */
        virtual void step2D(
            physics::RigidBodyAssembler& bodies,
            const Eigen::Vector2d& gravity,
            const double& time_step) const override;
    };

} // namespace time_stepper
} // namespace ccd
