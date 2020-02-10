#include "ccd/hash_grid.hpp"

#include <logger.hpp>

namespace ccd {

bool AABB::are_overlaping(const AABB& a, const AABB& b)
{
    // https://bit.ly/2ZP3tW4
    assert(a.dim == b.dim);
    return (abs(a.center.x() - b.center.x())
            <= (a.half_extent.x() + b.half_extent.x()))
        && (abs(a.center.y() - b.center.y())
            <= (a.half_extent.y() + b.half_extent.y()))
        && (a.dim == 2
            || abs(a.center.z() - b.center.z())
                <= (a.half_extent.z() + b.half_extent.z()));
};

void HashGrid::resize(Eigen::VectorXd min, Eigen::VectorXd max, double cellSize)
{
    clear();
    m_cellSize = cellSize;
    m_domainMin = min;
    m_domainMax = max;
    m_gridSize = int(std::ceil((max - min).maxCoeff() / m_cellSize));
}

/// @brief Compute an AABB around a given 2D mesh.
void calculate_mesh_extents(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    Eigen::VectorXd& lower_bound,
    Eigen::VectorXd& upper_bound)
{
    int dim = vertices.cols();
    Eigen::MatrixXd points(vertices.rows() + displacements.rows(), dim);
    points.topRows(vertices.rows()) = vertices;
    points.bottomRows(displacements.rows()) = vertices + displacements;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

/// @brief Compute the average edge length of a mesh.
double average_edge_length(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E)
{
    double avg = 0;
    for (unsigned i = 0; i < E.rows(); ++i) {
        avg += (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
        avg += ((V.row(E(i, 0)) + U.row(E(i, 0)))
                - (V.row(E(i, 1)) + U.row(E(i, 1))))
                   .norm();
    }
    return avg / (2 * E.rows());
}

/// @brief Compute the average displacement length.
double average_displacement_length(const Eigen::MatrixXd& displacements)
{
    return displacements.rowwise().norm().sum() / displacements.rows();
}

void HashGrid::resize(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    Eigen::VectorXd mesh_min, mesh_max;
    calculate_mesh_extents(vertices, displacements, mesh_min, mesh_max);
    this->resize(
        mesh_min.array() - inflation_radius,
        mesh_max.array() + inflation_radius,
        std::max(
            average_edge_length(vertices, displacements, edges),
            average_displacement_length(displacements))
                * 2
            + inflation_radius);
}

/// @brief Compute a AABB for a vertex moving through time (i.e. temporal edge).
void calculate_vertex_extents(
    const Eigen::VectorXd& v,
    const Eigen::VectorXd& u,
    Eigen::VectorXd& lower_bound,
    Eigen::VectorXd& upper_bound)
{
    Eigen::MatrixXd points(2, v.size());
    points.row(0) = v;
    points.row(1) = v + u;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addVertex(
    const Eigen::VectorXd& v,
    const Eigen::VectorXd& u,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorXd lower_bound, upper_bound;
    calculate_vertex_extents(v, u, lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_vertexItems); // Vertices have a negative id
}

void HashGrid::addVertices(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const double inflation_radius)
{
    tbb::parallel_for(0l, vertices.rows(), 1l, [&](long i) {
        addVertex(vertices.row(i), displacements.row(i), i, inflation_radius);
    });
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_edge_extents(
    const Eigen::VectorXd& vi,
    const Eigen::VectorXd& vj,
    const Eigen::VectorXd& ui,
    const Eigen::VectorXd& uj,
    Eigen::VectorXd& lower_bound,
    Eigen::VectorXd& upper_bound)
{
    Eigen::MatrixXd points(4, vi.size());
    points.row(0) = vi;
    points.row(1) = vj;
    points.row(2) = vi + ui;
    points.row(3) = vj + uj;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addEdge(
    const Eigen::VectorXd& vi,
    const Eigen::VectorXd& vj,
    const Eigen::VectorXd& ui,
    const Eigen::VectorXd& uj,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorXd lower_bound, upper_bound;
    calculate_edge_extents(vi, vj, ui, uj, lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_edgeItems); // Edges have a positive id
}

void HashGrid::addEdges(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    tbb::parallel_for(0l, edges.rows(), 1l, [&](long i) {
        addEdge(
            vertices.row(edges(i, 0)), vertices.row(edges(i, 1)),
            displacements.row(edges(i, 0)), displacements.row(edges(i, 1)), i,
            inflation_radius);
    });
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_face_extents(
    const Eigen::VectorXd& vi,
    const Eigen::VectorXd& vj,
    const Eigen::VectorXd& vk,
    const Eigen::VectorXd& ui,
    const Eigen::VectorXd& uj,
    const Eigen::VectorXd& uk,
    Eigen::VectorXd& lower_bound,
    Eigen::VectorXd& upper_bound)
{
    Eigen::MatrixXd points(6, vi.size());
    points.row(0) = vi;
    points.row(1) = vj;
    points.row(2) = vk;
    points.row(3) = vi + ui;
    points.row(4) = vj + uj;
    points.row(5) = vk + uk;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addFace(
    const Eigen::VectorXd& vi,
    const Eigen::VectorXd& vj,
    const Eigen::VectorXd& vk,
    const Eigen::VectorXd& ui,
    const Eigen::VectorXd& uj,
    const Eigen::VectorXd& uk,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorXd lower_bound, upper_bound;
    calculate_face_extents(vi, vj, vk, ui, uj, uk, lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_faceItems); // Faces have a positive id
}

void HashGrid::addFaces(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    tbb::parallel_for(0l, faces.rows(), 1l, [&](long i) {
        addFace(
            vertices.row(faces(i, 0)), vertices.row(faces(i, 1)),
            vertices.row(faces(i, 2)), displacements.row(faces(i, 0)),
            displacements.row(faces(i, 1)), displacements.row(faces(i, 2)), i,
            inflation_radius);
    });
}

void HashGrid::addElement(const AABB& aabb, const int id, HashItems& items)
{
    Eigen::VectorXi int_min =
        ((aabb.getMin() - m_domainMin) / m_cellSize).cast<int>();
    Eigen::VectorXi int_max =
        ((aabb.getMax() - m_domainMin) / m_cellSize).cast<int>();
    assert((int_min.array() <= int_max.array()).all());

    int min_z = int_min.size() == 3 ? int_min.z() : 0;
    int max_z = int_max.size() == 3 ? int_max.z() : 0;
    for (int x = int_min.x(); x <= int_max.x(); ++x) {
        for (int y = int_min.y(); y <= int_max.y(); ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                items.emplace_back(hash(x, y, z), id, aabb);
            }
        }
    }
}

template <typename T>
void getPairs(
    std::function<bool(int, int)> is_endpoint,
    std::function<bool(int, int)> is_same_group,
    HashItems& items0,
    HashItems& items1,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items0.begin(), items0.end());
    tbb::parallel_sort(items1.begin(), items1.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    int i = 0, j_start = 0;
    while (i < items0.size() && j_start < items1.size()) {
        const HashItem& item0 = items0[i];

        int j = j_start;
        while (j < items1.size()) {
            const HashItem& item1 = items1[j];

            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && !is_same_group(item0.id, item1.id)
                    && AABB::are_overlaping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break;
            }
            j++;
        }

        if (i == items0.size() - 1 || item0.key != items0[i + 1].key) {
            j_start = j;
        }
        i++;
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

template <typename T>
void getPairs(
    std::function<bool(int, int)> is_endpoint,
    std::function<bool(int, int)> is_same_group,
    HashItems& items,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items.begin(), items.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    for (int i = 0; i < items.size(); i++) {
        const HashItem& item0 = items[i];
        for (int j = i + 1; j < items.size(); j++) {
            const HashItem& item1 = items[j];
            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && !is_same_group(item0.id, item1.id)
                    && AABB::are_overlaping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break; // This avoids a brute force comparison
            }
        }
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

void HashGrid::getVertexEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates)
{
    auto is_endpoint = [&](int ei, int vi) {
        return edges(ei, 0) == vi || edges(ei, 1) == vi;
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int ei, int vi) {
        return check_groups
            && (group_ids(vi) == group_ids(edges(ei, 0))
                || group_ids(vi) == group_ids(edges(ei, 1)));
    };

    getPairs(
        is_endpoint, is_same_group, m_edgeItems, m_vertexItems, ev_candidates);
}

void HashGrid::getEdgeEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeEdgeCandidates& ee_candidates)
{
    auto is_endpoint = [&](int ei, int ej) {
        return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
            || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int ei, int ej) {
        return check_groups
            && (group_ids(edges(ei, 0)) == group_ids(edges(ej, 0))
                || group_ids(edges(ei, 0)) == group_ids(edges(ej, 1))
                || group_ids(edges(ei, 1)) == group_ids(edges(ej, 0))
                || group_ids(edges(ei, 1)) == group_ids(edges(ej, 1)));
    };

    getPairs(is_endpoint, is_same_group, m_edgeItems, ee_candidates);
}

void HashGrid::getFaceVertexPairs(
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    FaceVertexCandidates& fv_candidates)
{
    auto is_endpoint = [&](int fi, int vi) {
        return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi, 2);
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int fi, int vi) {
        return check_groups
            && (group_ids(vi) == group_ids(faces(fi, 0))
                || group_ids(vi) == group_ids(faces(fi, 1))
                || group_ids(vi) == group_ids(faces(fi, 2)));
    };

    getPairs(
        is_endpoint, is_same_group, m_faceItems, m_vertexItems, fv_candidates);
}

} // namespace ccd
