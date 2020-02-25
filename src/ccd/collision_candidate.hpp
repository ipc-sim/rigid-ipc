#pragma once

#include <string>
// #include <tbb/concurrent_vector.h>
#include <vector>

namespace ccd {

struct EdgeVertexCandidate {
    EdgeVertexCandidate(long edge_index, long vertex_index);

    bool operator==(const EdgeVertexCandidate& other) const;

    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    long edge_index;
    long vertex_index;
};

struct EdgeEdgeCandidate {
    EdgeEdgeCandidate(long edge0_index, long edge1_index);

    bool operator==(const EdgeEdgeCandidate& other) const;

    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const EdgeEdgeCandidate& other) const;

    long edge0_index;
    long edge1_index;
};

struct FaceVertexCandidate {
    FaceVertexCandidate(long face_index, long vertex_index);

    bool operator==(const FaceVertexCandidate& other) const;

    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    long face_index;
    long vertex_index;
};

struct Candidates {
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;

    size_t size() const
    {
        return ev_candidates.size() + ee_candidates.size()
            + fv_candidates.size();
    }

    void clear()
    {
        ev_candidates.clear();
        ee_candidates.clear();
        fv_candidates.clear();
    }
};

// struct ConcurrentCandidates {
//     tbb::concurrent_vector<EdgeVertexCandidate> ev_candidates;
//     tbb::concurrent_vector<EdgeEdgeCandidate> ee_candidates;
//     tbb::concurrent_vector<FaceVertexCandidate> fv_candidates;
//
//     size_t size() const
//     {
//         return ev_candidates.size() + ee_candidates.size()
//             + fv_candidates.size();
//     }
//
//     void clear()
//     {
//         ev_candidates.clear();
//         ee_candidates.clear();
//         fv_candidates.clear();
//     }
// };

} // namespace ccd
