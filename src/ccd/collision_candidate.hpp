#pragma once

#include <string>
#include <vector>

namespace ccd {

class Candidate {
public:
    virtual ~Candidate() = default;

    virtual bool operator==(const Candidate& other) const;

    /// @brief Compare EdgeVertexCandidates for sorting.
    virtual bool operator<(const Candidate& other) const;

    virtual std::string string() const;

protected:
    virtual inline int get_index0() const = 0;
    virtual inline int get_index1() const = 0;
};

class EdgeVertexCandidate : public Candidate {
public:
    EdgeVertexCandidate(long edge_index, long vertex_index);

    long edge_index;
    long vertex_index;

protected:
    virtual inline int get_index0() const override { return edge_index; };
    virtual inline int get_index1() const override { return vertex_index; };
};

class EdgeEdgeCandidate : public Candidate {
public:
    EdgeEdgeCandidate(long edge0_index, long edge1_index);

    virtual bool operator==(const EdgeEdgeCandidate& other) const;

    /// @brief Compare EdgeVertexCandidates for sorting.
    virtual bool operator<(const EdgeEdgeCandidate& other) const;

    long edge0_index;
    long edge1_index;

protected:
    virtual inline int get_index0() const override { return edge0_index; };
    virtual inline int get_index1() const override { return edge1_index; };
};

class FaceVertexCandidate : public Candidate {
public:
    FaceVertexCandidate(long face_index, long vertex_index);

    long face_index;
    long vertex_index;

protected:
    virtual inline int get_index0() const override { return face_index; };
    virtual inline int get_index1() const override { return vertex_index; };
};

typedef std::vector<EdgeVertexCandidate> EdgeVertexCandidates;
typedef std::vector<EdgeEdgeCandidate> EdgeEdgeCandidates;
typedef std::vector<FaceVertexCandidate> FaceVertexCandidates;

} // namespace ccd
