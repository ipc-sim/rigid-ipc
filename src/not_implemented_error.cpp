#include <FixingCollisions/not_implemented_error.hpp>

namespace ccd {
NotImplementedError::NotImplementedError(const char * error_message)
    : std::logic_error(error_message)
{
}

}
