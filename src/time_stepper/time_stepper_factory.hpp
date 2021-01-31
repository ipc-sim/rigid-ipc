#pragma once

#include <memory> // shared_ptr

#include <time_stepper/time_stepper.hpp>

namespace ipc::rigid {

class TimeStepperFactory {
public:
    static const TimeStepperFactory& factory();

    std::shared_ptr<TimeStepper>
    get_time_stepper(const std::string& name) const;

    std::shared_ptr<TimeStepper>
    get_default_time_stepper(int dim) const;

private:
    TimeStepperFactory();

    std::map<std::string, std::shared_ptr<TimeStepper>>
        time_steppers;
};

} // namespace ipc::rigid
