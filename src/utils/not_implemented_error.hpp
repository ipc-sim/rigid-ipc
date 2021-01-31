#pragma once

#include <stdexcept>
#include <string>

namespace ipc::rigid {

/// Error class for indicating features have not been implmented yet.
class NotImplementedError : public std::logic_error {
public:
    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    NotImplementedError(
        const char* error_message = "Functionality not yet implemented!")
        : std::logic_error(error_message) {};

    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    NotImplementedError(
        const std::string error_message = "Functionality not yet implemented!")
        : std::logic_error(error_message) {};
};

class DeprecatedError : public std::logic_error {
public:
    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    DeprecatedError(
        const char* error_message =
            "Functionality deprecated. Do not use this function!")
        : std::logic_error(error_message) {};

    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    DeprecatedError(
        const std::string error_message =
            "Functionality deprecated. Do not use this function!")
        : std::logic_error(error_message) {};
};

} // namespace ipc::rigid
