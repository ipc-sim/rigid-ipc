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
    m_cellSize = cellSize * 2.0;
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
    return avg / E.rows();
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
        (average_edge_length(vertices, displacements, edges)
         + average_displacement_length(displacements))
                / 2.0
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
        -(index + 1), m_vertexItems); // Vertices have a negative id
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
        index + 1, m_edgeItems); // Edges have a positive id
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
    points.row(1) = vk;
    points.row(2) = vi + ui;
    points.row(3) = vj + uj;
    points.row(3) = vk + uk;

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
        index + 1, m_faceItems); // Faces have a positive id
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
    Eigen::VectorXi int_min
        = ((aabb.getMin() - m_domainMin) / m_cellSize).cast<int>();
    Eigen::VectorXi int_max
        = ((aabb.getMax() - m_domainMin) / m_cellSize)
              .unaryExpr([](const double x) { return std::ceil(x); })
              .cast<int>();

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

void HashGrid::getVertexEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates)
{
    std::vector<HashItem> edge_items;
    std::vector<HashItem> vertex_items;

    bool check_groups = group_ids.size() > 0;

    // Combine the edge and vertex hash items.
    std::vector<HashItem> items;
    items.reserve(m_vertexItems.size() + m_edgeItems.size());
    items.insert(items.end(), m_vertexItems.begin(), m_vertexItems.end());
    items.insert(items.end(), m_edgeItems.begin(), m_edgeItems.end());

    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items.begin(), items.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    for (unsigned i = 0; i < items.size(); ++i) {
        HashItem& item = items[i];
        const int currH = item.key;
        const int currId = item.id;

        // read this element
        if (currId > 0) { // Edge elements have positive id
            item.id = currId - 1;
            edge_items.push_back(item);
        } else if (currId < 0) { // Vertex elements have negative id
            item.id = -currId - 1;
            vertex_items.push_back(item);
        } else {
            throw "Invalid id given to Hash!";
        }

        // CLOSE BUCKET:
        // if this is the last element, or next is from another bucket
        if ((i == items.size() - 1) || currH != items[i + 1].key) {
            // We are closing the bucket (key entry), so tally up all
            // vertex-edge pairs encountered in the bucket that just ended
            for (const HashItem& edge_item : edge_items) {
                const int& edge_id = edge_item.id;
                const AABB& edge_aabb = edge_item.aabb;

                for (const HashItem& vertex_item : vertex_items) {
                    const int& vertex_id = vertex_item.id;
                    const AABB& vertex_aabb = vertex_item.aabb;

                    bool is_endpoint = edges(edge_id, 0) == vertex_id
                        || edges(edge_id, 1) == vertex_id;
                    bool same_group = false;
                    if (check_groups) {
                        same_group = group_ids(vertex_id)
                            == group_ids(edges(edge_id, 0));
                    }
                    if (!is_endpoint && !same_group
                        && AABB::are_overlaping(edge_aabb, vertex_aabb)) {
                        ev_candidates.emplace_back(edge_id, vertex_id);
                    }
                }
            }

            edge_items.clear();
            vertex_items.clear();
        }
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(ev_candidates.begin(), ev_candidates.end());
    auto new_end = std::unique(ev_candidates.begin(), ev_candidates.end());
    ev_candidates.erase(new_end, ev_candidates.end());
}

void HashGrid::getEdgeEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    EdgeEdgeCandidates& ee_candidates)
{
    std::vector<HashItem> edge_items;

    bool check_groups = group_ids.size() > 0;

    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(m_edgeItems.begin(), m_edgeItems.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    for (unsigned i = 0; i < m_edgeItems.size(); ++i) {
        HashItem& item = m_edgeItems[i];
        const int currH = item.key;
        const int currId = item.id;

        // read this element
        if (currId > 0) { // Edge elements have positive id
            item.id = currId - 1;
            edge_items.push_back(item);
        } else {
            throw "Invalid edge id given to Hash!";
        }

        // CLOSE BUCKET:
        // if this is the last element, or next is from another bucket
        if ((i == m_edgeItems.size() - 1) || currH != m_edgeItems[i + 1].key) {
            // We are closing the bucket (key entry), so tally up all
            // edge-edge pairs encountered in the bucket that just ended
            for (int ei = 0; ei < edge_items.size(); ei++) {
                const HashItem& ei_item = edge_items[ei];
                const int& ei_id = ei_item.id;
                const AABB& ei_aabb = ei_item.aabb;

                for (int ej = ei; ej < edge_items.size(); ej++) {
                    const HashItem& ej_item = edge_items[ej];
                    const int& ej_id = ej_item.id;
                    const AABB& ej_aabb = ej_item.aabb;

                    bool has_common_endpoint
                        = edges(ei_id, 0) == edges(ej_id, 0)
                        || edges(ei_id, 0) == edges(ej_id, 1)
                        || edges(ei_id, 1) == edges(ej_id, 0)
                        || edges(ei_id, 1) == edges(ej_id, 1);
                    bool same_group = false;
                    if (check_groups) {
                        same_group = group_ids(edges(ei_id, 0))
                            == group_ids(edges(ej_id, 0));
                    }
                    if (!has_common_endpoint && !same_group
                        && AABB::are_overlaping(ei_aabb, ej_aabb)) {
                        ee_candidates.emplace_back(ei_id, ej_id);
                    }
                }
            }

            edge_items.clear();
        }
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(ee_candidates.begin(), ee_candidates.end());
    auto new_end = std::unique(ee_candidates.begin(), ee_candidates.end());
    ee_candidates.erase(new_end, ee_candidates.end());
}

void HashGrid::getFaceVertexPairs(
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    FaceVertexCandidates& fv_candidates)
{
    std::vector<HashItem> face_items;
    std::vector<HashItem> vertex_items;

    bool check_groups = group_ids.size() > 0;

    // Combine the face and vertex hash items.
    std::vector<HashItem> items;
    items.reserve(m_vertexItems.size() + m_faceItems.size());
    items.insert(items.end(), m_vertexItems.begin(), m_vertexItems.end());
    items.insert(items.end(), m_faceItems.begin(), m_faceItems.end());

    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items.begin(), items.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-face pairs with the same key
    for (unsigned i = 0; i < items.size(); ++i) {
        HashItem& item = items[i];
        const int currH = item.key;
        const int currId = item.id;

        // read this element
        if (currId > 0) { // Face elements have positive id
            item.id = currId - 1;
            face_items.push_back(item);
        } else if (currId < 0) { // Vertex elements have negative id
            item.id = -currId - 1;
            vertex_items.push_back(item);
        } else {
            throw "Invalid id given to Hash!";
        }

        // CLOSE BUCKET:
        // if this is the last element, or next is from another bucket
        if ((i == items.size() - 1) || currH != items[i + 1].key) {
            // We are closing the bucket (key entry), so tally up all
            // vertex-face pairs encountered in the bucket that just ended
            for (const HashItem& face_item : face_items) {
                const int& face_id = face_item.id;
                const AABB& face_aabb = face_item.aabb;

                for (const HashItem& vertex_item : vertex_items) {
                    const int& vertex_id = vertex_item.id;
                    const AABB& vertex_aabb = vertex_item.aabb;

                    bool is_endpoint = faces(face_id, 0) == vertex_id
                        || faces(face_id, 1) == vertex_id
                        || faces(face_id, 2) == vertex_id;
                    bool same_group = false;
                    if (check_groups) {
                        same_group = group_ids(vertex_id)
                            == group_ids(faces(face_id, 0));
                    }
                    if (!is_endpoint && !same_group
                        && AABB::are_overlaping(face_aabb, vertex_aabb)) {
                        fv_candidates.emplace_back(face_id, vertex_id);
                    }
                }
            }

            face_items.clear();
            vertex_items.clear();
        }
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(fv_candidates.begin(), fv_candidates.end());
    auto new_end = std::unique(fv_candidates.begin(), fv_candidates.end());
    fv_candidates.erase(new_end, fv_candidates.end());
}

} // namespace ccd
