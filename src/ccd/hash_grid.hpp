#pragma once

#include <iomanip>
#include <list>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
namespace ccd {

struct EdgeVertexCandidate {
    long edge_index;
    long vertex_index;

    EdgeVertexCandidate(long edge_index, long vertex_index)
        : edge_index(edge_index)
        , vertex_index(vertex_index)
    {
    }

    bool operator==(const EdgeVertexCandidate& other) const
    {
        return this->edge_index == other.edge_index
            && this->vertex_index == other.vertex_index;
    }
};

typedef std::vector<EdgeVertexCandidate> EdgeVertexCandidates;
typedef std::unordered_set<EdgeVertexCandidate> EdgeVertexCandidateSet;

} // namespace ccd

namespace std {
template <> struct hash<ccd::EdgeVertexCandidate> {
    inline size_t operator()(const ccd::EdgeVertexCandidate& ev_candidate) const
    {
        // std::hash<size_t> int_hasher;
        // return int_hasher(v.first) ^ int_hasher(v.second);
        return ev_candidate.edge_index * 31 + ev_candidate.vertex_index;
    }
};
} // namespace std

namespace ccd {

/// @brief An entry into the hash grid as a (key, value) pair.
class HashItem {
public:
    int key; /// @brief The key of the item.
    int id;  /// @brief The value of the item.

    /// @breif Construct a hash item as a (key, value) pair.
    HashItem(int key, int id)
        : key(key)
        , id(id)
    {
    }

    /// @brief Compare HashItems by their keys for sorting.
    bool operator<(const HashItem& other) const { return key < other.key; }
};

/// @brief Axis aligned bounding-box of some type
template <class T> class AABBT {
public:
    T min;
    T max;

    AABBT()
        : AABBT(T(0, 0), T(0, 0))
    {
    }

    AABBT(const T& min, const T& max)
        : min(min)
        , max(max)
    {
    }

    virtual ~AABBT() {}
};

typedef AABBT<Eigen::Vector2i> AABBi;
typedef AABBT<Eigen::Vector2d> AABBd;

class HashGrid {
public:
    void resize(Eigen::Vector2d mn, Eigen::Vector2d mx, double cellSize);
    void resize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i edges);

    void addVertex(
        const Eigen::Vector2d& v, const Eigen::Vector2d& u, const int index);

    void addVertices(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements);

    void addEdge(const Eigen::Vector2d& vi, const Eigen::Vector2d& vj,
        const Eigen::Vector2d& ui, const Eigen::Vector2d& uj, const int index);

    void addEdges(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i edges);

    void getVertexEdgePairs(
        const Eigen::MatrixX2i& edges, EdgeVertexCandidates& ev_candidates);

protected:
    void addElement(Eigen::Vector2d xmin, Eigen::Vector2d xmax, int id);

    AABBi makeAABBi(Eigen::Vector2d mn, Eigen::Vector2d mx);

    void sort();

    void add(AABBi aabbi, int id);
    void add(Eigen::Vector2i p, int id);

    int hash(Eigen::Vector2i p);
    void clear();

    HashItem& get(unsigned int i);

protected:
    double m_cellSize;
    int m_gridSize;
    Eigen::Vector2d m_domainMin;
    Eigen::Vector2d m_domainMax;

    std::vector<HashItem> m_hash;
};

} // namespace ccd
