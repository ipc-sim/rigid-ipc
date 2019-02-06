#include <FixingCollisions/impact.hpp>

namespace ccd {

Impact::Impact(int vertex_index, int edge_index, double time, double alpha)
    : vertex_index(vertex_index)
    , edge_index(edge_index)
    , time(time)
    , alpha(alpha)
{
}

Impact::~Impact() {}

}
