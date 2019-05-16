// Hash.h
//
// This is a slightly modified version of code written by
// Daniele Panozzo
//

#ifndef _HASH_H_
#define _HASH_H_

#include <iomanip>
#include <list>
#include <vector>

#include <Eigen/Core>

typedef std::pair<size_t, size_t> Candidate;
typedef std::vector<Candidate> Candidates;
typedef std::vector<Candidate>::iterator CandidatesIterator;

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

    /// @brief Construct a hash item as a empty pair.
    HashItem() {}

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

    void addElement(Eigen::Vector2d xmin, Eigen::Vector2d xmax, int id);

    void getVertexEdgePairs(const Eigen::MatrixX2i& edges, Candidates& hits);

protected:
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

#endif
