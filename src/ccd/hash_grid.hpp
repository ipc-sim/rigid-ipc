#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <list>
#include <tbb/parallel_sort.h>
#include <unordered_set>
#include <vector>

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
        // https://www.techiedelight.com/use-std-pair-key-std-unordered_map-cpp/
        std::hash<size_t> int_hasher;
        return int_hasher(ev_candidate.edge_index)
            ^ int_hasher(ev_candidate.vertex_index);
    }
};
} // namespace std

namespace ccd {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    AABB()
        : AABB(Eigen::Vector2d(0, 0), Eigen::Vector2d(0, 0))
    {
    }

    AABB(const Eigen::Vector2d& min, const Eigen::Vector2d& max)
        : min(min)
        , max(max)
    {
        half_width = (max.x() - min.x()) / 2.0;
        half_height = (max.y() - min.y()) / 2.0;
        center = min + Eigen::Vector2d(half_width, half_height);
    }

    virtual ~AABB() {}

    static bool are_overlaping(const AABB& a, const AABB& b);

    Eigen::Vector2d getMin() const { return min; }
    Eigen::Vector2d getMax() const { return max; }
    double getHalfWidth() const { return half_width; }
    double getHalfHeight() const { return half_height; }
    Eigen::Vector2d getCenter() const { return center; }

private:
    Eigen::Vector2d min;
    Eigen::Vector2d max;
    double half_width, half_height;
    Eigen::Vector2d center;
};

/// @brief An entry into the hash grid as a (key, value) pair.
class HashItem {
public:
    int key;   /// @brief The key of the item.
    int id;    /// @brief The value of the item.
    AABB aabb; /// @brief The axis-aligned bounding box of the element

    /// @breif Construct a hash item as a (key, value) pair.
    HashItem(int key, int id, const AABB aabb)
        : key(key)
        , id(id)
        , aabb(aabb)
    {
    }

    /// @brief Compare HashItems by their keys for sorting.
    bool operator<(const HashItem& other) const { return key < other.key; }
};

class HashGrid {
public:
    void resize(Eigen::Vector2d mn, Eigen::Vector2d mx, double cellSize);
    void resize(
        const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements,
        const Eigen::MatrixX2i edges,
        const double inflation_radius = 0.0);

    /// @brief Add a vertex as a AABB containing the time swept edge.
    void addVertex(
        const Eigen::Vector2d& v,
        const Eigen::Vector2d& u,
        const int index,
        const double inflation_radius = 0.0);

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVertices(
        const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addEdge(
        const Eigen::Vector2d& vi,
        const Eigen::Vector2d& vj,
        const Eigen::Vector2d& ui,
        const Eigen::Vector2d& uj,
        const int index,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addEdges(
        const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements,
        const Eigen::MatrixX2i& edges,
        const double inflation_radius = 0.0);

    /// @brief Compute the candidate edge-vertex intersections.
    void getVertexEdgePairs(
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXi& group_ids,
        EdgeVertexCandidates& ev_candidates);

protected:
    /// @brief Add an AABB of the extents to the hash grid.
    void addElement(const AABB& aabb, const int id);

    /// @brief Sort all hash items.
    inline void sort() { tbb::parallel_sort(m_hash.begin(), m_hash.end()); }

    /// @brief Create the hash of a cell location.
    inline int hash(const int& x, const int& y) const
    {
        return y * m_gridSize + x;
    }

    /// @brief Clear the hash grid.
    inline void clear() { m_hash.clear(); }

    /// @brief Get the item with the given id.
    inline HashItem& get(unsigned int i) { return m_hash[i]; }

protected:
    double m_cellSize;
    int m_gridSize;
    Eigen::Vector2d m_domainMin;
    Eigen::Vector2d m_domainMax;

    std::vector<HashItem> m_hash;
};

} // namespace ccd
