#include "collision_candidate.hpp"

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

EdgeVertexCandidate::EdgeVertexCandidate(long edge_index, long vertex_index)
    : edge_index(edge_index)
    , vertex_index(vertex_index)
{
}

EdgeEdgeCandidate::EdgeEdgeCandidate(long edge0_index, long edge1_index)
    : edge0_index(edge0_index)
    , edge1_index(edge1_index)
{
    // This should never happen
    assert(edge0_index != edge1_index);
}

FaceVertexCandidate::FaceVertexCandidate(long face_index, long vertex_index)
    : face_index(face_index)
    , vertex_index(vertex_index)
{
}

} // namespace ccd
