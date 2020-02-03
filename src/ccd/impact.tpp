#pragma once

#include <ccd/impact.hpp>

namespace ccd {

// Compare two edge-vertex impacts to determine if impact0 comes before impact1.
template <typename Impact>
bool compare_impacts_by_time(Impact impact1, Impact impact2)
{
    return impact1.time < impact2.time;
}

} // namespace ccd
