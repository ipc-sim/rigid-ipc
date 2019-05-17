#pragma once

#include <iomanip>
#include <list>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
namespace std {
template <> struct hash<std::pair<size_t, size_t>> {
    inline size_t operator()(const std::pair<int, int>& v) const
    {
        // std::hash<size_t> int_hasher;
        // return int_hasher(v.first) ^ int_hasher(v.second);
        return v.first * 31 + v.second;
    }
};
} // namespace std

namespace ccd {

typedef std::pair<size_t, size_t> Candidate;
typedef std::unordered_set<Candidate> Candidates;
// typedef std::unordered_set<Candidate, boost::hash<Candidate>>::iterator
// CandidatesIterator;

/// @brief An entry into the hash grid as a (key, value) pair.
class HashItem {
public:
    int key; /// @brief The key of the item.
    int id;  /// @brief The value of the item.

    /// @breif Construct a hash item as a (key, value) pair.
    HashItem(int k, int i)
    {
        key = k;
        id = i;
    }

    /// @brief Compare HashItems by their keys for sorting.
    bool operator<(const HashItem& hi) const { return key < hi.key; }
};

/// @brief Axis aligned bounding-box of some type
template <class T> class AABBT {
public:
    AABBT(const T& min, const T& max)
    {
        m_min = min;
        m_max = max;
    }

    AABBT()
    {
        m_min = T(0, 0);
        m_max = T(0, 0);
    }

    ~AABBT() {}

    T& getMin() { return m_min; }

    T& getMax() { return m_max; }

    T m_min;
    T m_max;
};

typedef Eigen::Matrix<int, 2, 1> Point2i;
typedef Eigen::Matrix<double, 2, 1> Point2d;
typedef AABBT<Point2i> AABBi;
typedef AABBT<Point2d> AABBd;

class Hash {
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

    void getVertexEdgePairs(const Eigen::MatrixX2i& edges, Candidates& hits);

protected:
    void addElement(Eigen::Vector2d xmin, Eigen::Vector2d xmax, int id);

    AABBi makeAABBi(Eigen::Vector2d mn, Eigen::Vector2d mx);

    void sort();

    void add(AABBi aabbi, int id);
    void add(Point2i p, int id);

    int hash(Point2i p);
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
