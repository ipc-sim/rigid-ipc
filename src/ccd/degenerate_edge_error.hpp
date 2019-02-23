#ifndef DEGENERATE_EDGE_ERROR_H
#define DEGENERATE_EDGE_ERROR_H

#include <stdexcept>

namespace ccd {

/// Error class for indicating a edge is in fact just a single point.
class DegenerateEdgeError : public std::logic_error {
public:
    /**
     * Create a degenerate edge error.

     * @param error_message Message describing the error.
     */
    DegenerateEdgeError(
        const char* error_message = "Edge is a singular point!");
};

}

#endif // DEGENERATE_EDGE_ERROR_H
