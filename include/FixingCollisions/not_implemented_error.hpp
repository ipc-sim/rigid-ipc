#ifndef NOT_IMPLEMENTED_ERROR_H
#define NOT_IMPLEMENTED_ERROR_H

#include <stdexcept>
namespace ccd {
class NotImplementedError : public std::logic_error {
public:
    NotImplementedError(const char * error_message = "Functionality not yet implemented!");
};
}
#endif // NOT_IMPLEMENTED_ERROR_H
