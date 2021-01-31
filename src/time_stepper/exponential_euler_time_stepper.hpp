#pragma once

#include <time_stepper/time_stepper.hpp>

namespace ipc::rigid {

class ExponentialEulerTimeStepper : public TimeStepper {
public:
    virtual ~ExponentialEulerTimeStepper() = default;

    static std::string default_name() { return "exponential_euler"; }
    virtual std::string name() const override
    {
        return ExponentialEulerTimeStepper::default_name();
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
        RigidBody& body,
        const Eigen::Vector3d& gravity,
        const double& time_step) const override;
};

} // namespace ipc::rigid
