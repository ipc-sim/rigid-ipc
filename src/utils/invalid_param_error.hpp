#pragma once

#include <stdexcept>

namespace ccd {

/// Error class for indicating features have not been implmented yet.
class InvalidParameterError : public std::logic_error {
public:
    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    InvalidParameterError(
        const char* error_message = "InvalidParameter")
        : std::logic_error(error_message) {};
};

} // namespace ccd
