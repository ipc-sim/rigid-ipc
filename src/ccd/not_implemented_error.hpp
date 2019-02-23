#ifndef NOT_IMPLEMENTED_ERROR_H
#define NOT_IMPLEMENTED_ERROR_H

#include <stdexcept>

namespace ccd {

/// Error class for indicating features have not been implmented yet.
class NotImplementedError : public std::logic_error {
public:
    /**
     * Create a not implemented error.
     *
     * @param error_message Message describing what has not been implmented.
     */
    NotImplementedError(
        const char* error_message = "Functionality not yet implemented!");
};

}
#endif // NOT_IMPLEMENTED_ERROR_H
