#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <list>
#include <tbb/tbb.h>
#include <vector>

#include <ccd/collision_candidate.hpp>

namespace ccd {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    AABB() {}

    AABB(const Eigen::VectorXd& min, const Eigen::VectorXd& max)
        : min(min)
        , max(max)
    {
        half_extent = (max - min) / 2;
        center = min + half_extent;
        dim = min.size();
        assert(max.size() == dim);
    }

    virtual ~AABB() {}

    static bool are_overlaping(const AABB& a, const AABB& b);

    inline const Eigen::VectorXd& getMin() const { return min; }
    inline const Eigen::VectorXd& getMax() const { return max; }
    inline const Eigen::VectorXd& getHalfExtent() const { return half_extent; }
    inline const Eigen::VectorXd& getCenter() const { return center; }

private:
    Eigen::VectorXd min;
    Eigen::VectorXd max;
    Eigen::VectorXd half_extent;
    Eigen::VectorXd center;
    int dim;
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

typedef tbb::concurrent_vector<HashItem> HashItems;

class HashGrid {
public:
    void resize(Eigen::VectorXd min, Eigen::VectorXd max, double cellSize);

    void resize(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXd& displacements,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add a vertex as a AABB containing the time swept edge.
    void addVertex(
        const Eigen::VectorXd& v,
        const Eigen::VectorXd& u,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVertices(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXd& displacements,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addEdge(
        const Eigen::VectorXd& vi,
        const Eigen::VectorXd& vj,
        const Eigen::VectorXd& ui,
        const Eigen::VectorXd& uj,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addEdges(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXd& displacements,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addFace(
        const Eigen::VectorXd& vi,
        const Eigen::VectorXd& vj,
        const Eigen::VectorXd& vk,
        const Eigen::VectorXd& ui,
        const Eigen::VectorXd& uj,
        const Eigen::VectorXd& uk,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addFaces(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXd& displacements,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0.0);

    /// @brief Compute the candidate edge-vertex candidate collisisons.
    void getVertexEdgePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::VectorXi& group_ids,
        EdgeVertexCandidates& ev_candidates);

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getEdgeEdgePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::VectorXi& group_ids,
        EdgeEdgeCandidates& ee_candidates);

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getFaceVertexPairs(
        const Eigen::MatrixXi& faces,
        const Eigen::VectorXi& group_ids,
        FaceVertexCandidates& fv_candidates);

protected:
    /// @brief Add an AABB of the extents to the hash grid.
    void addElement(const AABB& aabb, const int id, HashItems& items);

    /// @brief Create the hash of a cell location.
    inline long hash(int x, int y, int z) const
    {
        return z + m_gridSize * (y * m_gridSize + x);
    }

    /// @brief Clear the hash grid.
    inline void clear()
    {
        m_vertexItems.clear();
        m_edgeItems.clear();
        m_faceItems.clear();
    }

protected:
    double m_cellSize;
    int m_gridSize;
    Eigen::VectorXd m_domainMin;
    Eigen::VectorXd m_domainMax;

    HashItems m_vertexItems;
    HashItems m_edgeItems;
    HashItems m_faceItems;
};

} // namespace ccd
