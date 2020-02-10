#include "collision_candidate.hpp"

#include <fmt/format.h>

namespace ccd {

bool Candidate::operator==(const Candidate& other) const
{
    return this->get_index0() == other.get_index0()
        && this->get_index1() == other.get_index1();
}

bool Candidate::operator<(const Candidate& other) const
{
    if (this->get_index0() == other.get_index0()) {
        return this->get_index1() < other.get_index1();
    }
    return this->get_index0() < other.get_index0();
}

std::string Candidate::string() const
{
    return fmt::format("({:d}, {:d})", this->get_index0(), this->get_index1());
}

EdgeVertexCandidate::EdgeVertexCandidate(long edge_index, long vertex_index)
    : edge_index(edge_index)
    , vertex_index(vertex_index)
{
}

EdgeEdgeCandidate::EdgeEdgeCandidate(long edge0_index, long edge1_index)
    : edge0_index(edge0_index)
    , edge1_index(edge1_index)
{
}

bool EdgeEdgeCandidate::operator==(const EdgeEdgeCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->edge0_index == other.edge0_index
            && this->edge1_index == other.edge1_index)
        || (this->edge0_index == other.edge1_index
            && this->edge1_index == other.edge0_index);
}

bool EdgeEdgeCandidate::operator<(const EdgeEdgeCandidate& other) const
{
    long this_min = std::min(this->edge0_index, this->edge1_index);
    long other_min = std::min(other.edge0_index, other.edge1_index);
    if (this_min == other_min) {
        return std::max(this->edge0_index, this->edge1_index)
            < std::max(other.edge0_index, other.edge1_index);
    }
    return this_min < other_min;
}

FaceVertexCandidate::FaceVertexCandidate(long face_index, long vertex_index)
    : face_index(face_index)
    , vertex_index(vertex_index)
{
}

} // namespace ccd
