#include "time_stepper_factory.hpp"

#include <time_stepper/dmv_time_stepper.hpp>
#include <time_stepper/exponential_euler_time_stepper.hpp>
#include <time_stepper/sympletic_euler_time_stepper.hpp>
#include <time_stepper/verlet_time_stepper.hpp>

namespace ipc::rigid {

const TimeStepperFactory& TimeStepperFactory::factory()
{
    static TimeStepperFactory instance;
    return instance;
}

TimeStepperFactory::TimeStepperFactory()
{
    // 2D
    time_steppers.emplace(
        SympleticEulerTimeStepper::default_name(),
        std::make_shared<SympleticEulerTimeStepper>());
    time_steppers.emplace(
        VerletTimeStepper::default_name(),
        std::make_shared<VerletTimeStepper>());
    // 3D
    time_steppers.emplace(
        ExponentialEulerTimeStepper::default_name(),
        std::make_shared<ExponentialEulerTimeStepper>());
    time_steppers.emplace(
        DMVTimeStepper::default_name(), std::make_shared<DMVTimeStepper>());
}

std::shared_ptr<TimeStepper>
TimeStepperFactory::get_time_stepper(const std::string& name) const
{
    auto it = time_steppers.find(name);
    if (it == time_steppers.end()) {
        spdlog::error("Invalid choice of time-stepper: {}", name);
    }
    assert(it != time_steppers.end());
    return it->second;
}

std::shared_ptr<TimeStepper>
TimeStepperFactory::get_default_time_stepper(int dim) const
{
    assert(dim == 2 || dim == 3);
    return dim == 2
        ? get_time_stepper(SympleticEulerTimeStepper::default_name())
        : get_time_stepper(DMVTimeStepper::default_name());
}

} // namespace ipc::rigid
