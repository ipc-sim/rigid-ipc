#include <ccd/hash_grid.hpp>
#include <logger.hpp>

namespace ccd {

void HashGrid::resize(Eigen::Vector2d mn, Eigen::Vector2d mx, double cellSize)
{
    clear();
    m_cellSize = cellSize * 2.0;
    m_domainMin = mn;
    m_domainMax = mx;
    m_gridSize = int(
        std::ceil(std::max(mx.x() - mn.x(), mx.y() - mn.y()) / m_cellSize));
}

/// @brief Compute an AABB around a given 2D mesh.
void calculate_mesh_extents(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    Eigen::Vector2d& lower_bound,
    Eigen::Vector2d& upper_bound)
{
    Eigen::MatrixXd points(vertices.rows() + displacements.rows(), 2);
    points.block(0, 0, vertices.rows(), vertices.cols()) = vertices;
    points.block(vertices.rows(), 0, displacements.rows(), displacements.cols())
        = vertices + displacements;

    lower_bound.x() = points.col(0).minCoeff();
    lower_bound.y() = points.col(1).minCoeff();
    upper_bound.x() = points.col(0).maxCoeff();
    upper_bound.y() = points.col(1).maxCoeff();
}

/// @brief Compute the average edge length of a mesh.
double average_edge_length(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const Eigen::MatrixX2i& edges)
{
    double sum = 0;
    for (int i = 0; i < edges.rows(); i++) {
        sum += (vertices.row(edges(i, 1)) - vertices.row(edges(i, 0))).norm();
        sum += ((vertices.row(edges(i, 1)) + displacements.row(edges(i, 1)))
            - (vertices.row(edges(i, 0)) + displacements.row(edges(i, 0))))
                   .norm();
    }
    return sum / (2 * edges.rows());
}

/// @brief Compute the average displacement length.
double average_displacement_length(const Eigen::MatrixX2d& displacements)
{
    double sum = 0;
    for (int i = 0; i < displacements.rows(); i++) {
        sum += displacements.row(i).norm();
    }
    return sum / displacements.rows();
}

void HashGrid::resize(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const Eigen::MatrixX2i edges,
    const double inflation_radius)
{
    static Eigen::Vector2d mesh_min, mesh_max;
    calculate_mesh_extents(vertices, displacements, mesh_min, mesh_max);
    this->resize(mesh_min.array() - inflation_radius,
        mesh_max.array() + inflation_radius,
        average_edge_length(vertices, displacements, edges) + inflation_radius
            + average_displacement_length(displacements));
}

/// @brief Compute a AABB for a vertex moving through time (i.e. temporal edge).
void calculate_vertex_extents(const Eigen::Vector2d& v,
    const Eigen::Vector2d& u,
    Eigen::Vector2d& lower_bound,
    Eigen::Vector2d& upper_bound)
{
    static Eigen::Matrix<double, 2, 2> points;
    points.row(0) = v;
    points.row(1) = v + u;

    lower_bound.x() = points.col(0).minCoeff();
    lower_bound.y() = points.col(1).minCoeff();
    upper_bound.x() = points.col(0).maxCoeff();
    upper_bound.y() = points.col(1).maxCoeff();
}

void HashGrid::addVertex(const Eigen::Vector2d& v,
    const Eigen::Vector2d& u,
    const int index,
    const double inflation_radius)
{
    static Eigen::Vector2d lower_bound, upper_bound;
    calculate_vertex_extents(v, u, lower_bound, upper_bound);
    this->addElement(lower_bound.array() - inflation_radius,
        upper_bound.array() + inflation_radius,
        -(index + 1)); // Vertices have a negative id
}

void HashGrid::addVertices(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const double inflation_radius)
{
    for (int i = 0; i < vertices.rows(); i++) {
        this->addVertex(
            vertices.row(i), displacements.row(i), i, inflation_radius);
    }
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_edge_extents(const Eigen::Vector2d& vi,
    const Eigen::Vector2d& vj,
    const Eigen::Vector2d& ui,
    const Eigen::Vector2d& uj,
    Eigen::Vector2d& lower_bound,
    Eigen::Vector2d& upper_bound)
{
    static Eigen::Matrix<double, 4, 2> points;
    points.row(0) = vi;
    points.row(1) = vj;
    points.row(2) = vi + ui;
    points.row(3) = vj + uj;

    lower_bound.x() = points.col(0).minCoeff();
    lower_bound.y() = points.col(1).minCoeff();
    upper_bound.x() = points.col(0).maxCoeff();
    upper_bound.y() = points.col(1).maxCoeff();
}

void HashGrid::addEdge(const Eigen::Vector2d& vi,
    const Eigen::Vector2d& vj,
    const Eigen::Vector2d& ui,
    const Eigen::Vector2d& uj,
    const int index,
    const double inflation_radius)
{
    static Eigen::Vector2d lower_bound, upper_bound;
    calculate_edge_extents(vi, vj, ui, uj, lower_bound, upper_bound);
    this->addElement(lower_bound.array() - inflation_radius,
        upper_bound.array() + inflation_radius,
        index + 1); // Edges have a positive id
}

void HashGrid::addEdges(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const Eigen::MatrixX2i& edges,
    const double inflation_radius)
{
    for (int i = 0; i < edges.rows(); i++) {
        this->addEdge(vertices.row(edges(i, 0)), vertices.row(edges(i, 1)),
            displacements.row(edges(i, 1)), displacements.row(edges(i, 1)), i,
            inflation_radius);
    }
}

void HashGrid::addElement(const Eigen::Vector2d& lower_bound,
    const Eigen::Vector2d& upper_bound,
    const int id)
{

    AABBi aabbi;
    aabbi.min = ((lower_bound - m_domainMin) / m_cellSize).cast<int>();
    aabbi.max = ((upper_bound - m_domainMin) / m_cellSize)
                    .unaryExpr([](const double x) { return std::ceil(x); })
                    .cast<int>();

    for (int x = aabbi.min.x(); x <= aabbi.max.x(); ++x) {
        for (int y = aabbi.min.y(); y <= aabbi.max.y(); ++y) {
            m_hash.push_back(HashItem(hash(x, y), id));
        }
    }
}

void HashGrid::getVertexEdgePairs(const Eigen::MatrixX2i& edges,
    const Eigen::VectorXi& group_ids,
    EdgeVertexCandidates& ev_candidates)
{
    EdgeVertexCandidateSet unique_ev_candidates;

    std::vector<int> edge_ids;
    std::vector<int> vertex_ids;

    bool check_groups = group_ids.size() > 0;

    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    this->sort();

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    for (unsigned i = 0; i < m_hash.size(); ++i) {
        const int& currH = get(i).key;
        const int& currId = get(i).id;

        // read this element
        if (currId > 0) { // Edge elements have positive id
            edge_ids.push_back(currId - 1);
        } else if (currId < 0) { // Vertex elements have negative id
            vertex_ids.push_back(-currId - 1);
        } else {
            throw "Invalid id given to Hash!";
        }

        // CLOSE BUCKET:
        // if this is the last element, or next is from another bucket
        if ((i == m_hash.size() - 1) || currH != get(i + 1).key) {
            // We are closing the bucket (key entry), so tally up all
            // vertex-edge pairs encountered in the bucket that just ended
            for (const int& edge_id : edge_ids) {
                for (const int& vertex_id : vertex_ids) {
                    bool is_endpoint = edges(edge_id, 0) == vertex_id
                        || edges(edge_id, 1) == vertex_id;
                    bool same_group = false;
                    if (check_groups) {
                        same_group = group_ids(vertex_id)
                            == group_ids(edges(edge_id, 0));
                    }
                    if (!is_endpoint && !same_group) {
                        unique_ev_candidates.insert(
                            EdgeVertexCandidate(edge_id, vertex_id));
                    }
                }
            }

            edge_ids.clear();
            vertex_ids.clear();
        }
    }

    // Copy the unique candidates over to the output vector of candidates
    std::copy(unique_ev_candidates.begin(), unique_ev_candidates.end(),
        std::back_inserter(ev_candidates));
}

} // namespace ccd
