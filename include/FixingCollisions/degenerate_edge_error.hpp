#ifndef DEGENERATE_EDGE_ERROR_H
#define DEGENERATE_EDGE_ERROR_H

#include <stdexcept>
namespace ccd {
class DegenerateEdgeError : public std::logic_error {
public:
    DegenerateEdgeError(
        const char* error_message = "Edge is a singular point!");
};
}
#endif // DEGENERATE_EDGE_ERROR_H
