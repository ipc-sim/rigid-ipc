#include <FixingCollisions/degenerate_edge_error.hpp>

namespace ccd {
DegenerateEdgeError::DegenerateEdgeError(const char* error_message)
    : std::logic_error(error_message)
{
}

}
