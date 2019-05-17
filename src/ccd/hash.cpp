#include <ccd/hash.hpp>

namespace ccd {

void Hash::resize(Eigen::Vector2d mn, Eigen::Vector2d mx, double cellSize)
{
    clear();
    m_cellSize = cellSize * 2.0;
    m_domainMin = mn;
    m_domainMax = mx;
    m_gridSize = int(
        std::ceil(std::max(mx.x() - mn.x(), mx.y() - mn.y()) / m_cellSize));
}

void calculate_mesh_extents(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, Eigen::Vector2d& AABB_min,
    Eigen::Vector2d& AABB_max)
{
    Eigen::MatrixXd points(vertices.rows() + displacements.rows(), 2);
    points.block(0, 0, vertices.rows(), vertices.cols()) = vertices;
    points.block(vertices.rows(), 0, displacements.rows(), displacements.cols())
        = vertices + displacements;

    AABB_min.x() = points.col(0).minCoeff();
    AABB_min.y() = points.col(1).minCoeff();
    AABB_max.x() = points.col(0).maxCoeff();
    AABB_max.y() = points.col(1).maxCoeff();
}

double average_edge_length(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges)
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

double average_displacement_length(const Eigen::MatrixX2d& displacements)
{
    double sum = 0;
    for (int i = 0; i < displacements.rows(); i++) {
        sum += displacements.row(i).norm();
    }
    return sum / displacements.rows();
}

void Hash::resize(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i edges)
{
    static Eigen::Vector2d mesh_min, mesh_max;
    calculate_mesh_extents(vertices, displacements, mesh_min, mesh_max);
    this->resize(mesh_min, mesh_max,
        average_edge_length(vertices, displacements, edges)
            + average_displacement_length(displacements));
}

void calculate_vertex_extents(const Eigen::Vector2d& v,
    const Eigen::Vector2d& u, Eigen::Vector2d& AABB_min,
    Eigen::Vector2d& AABB_max)
{
    static Eigen::Matrix<double, 2, 2> points;
    points.row(0) = v;
    points.row(1) = v + u;

    AABB_min.x() = points.col(0).minCoeff();
    AABB_min.y() = points.col(1).minCoeff();
    AABB_max.x() = points.col(0).maxCoeff();
    AABB_max.y() = points.col(1).maxCoeff();
}

void Hash::addVertex(
    const Eigen::Vector2d& v, const Eigen::Vector2d& u, const int index)
{
    static Eigen::Vector2d vertex_min, vertex_max;
    calculate_vertex_extents(v, u, vertex_min, vertex_max);
    this->addElement(
        vertex_min, vertex_max, -(index + 1)); // Vertices have a negative id
}

void Hash::addVertices(
    const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2d& displacements)
{
    for (int i = 0; i < vertices.rows(); i++) {
        this->addVertex(vertices.row(i), displacements.row(i), i);
    }
}

void calculate_edge_extents(const Eigen::Vector2d& vi,
    const Eigen::Vector2d& vj, const Eigen::Vector2d& ui,
    const Eigen::Vector2d& uj, Eigen::Vector2d& AABB_min,
    Eigen::Vector2d& AABB_max)
{
    static Eigen::Matrix<double, 4, 2> points;
    points.row(0) = vi;
    points.row(1) = vj;
    points.row(2) = vi + ui;
    points.row(3) = vj + uj;

    AABB_min.x() = points.col(0).minCoeff();
    AABB_min.y() = points.col(1).minCoeff();
    AABB_max.x() = points.col(0).maxCoeff();
    AABB_max.y() = points.col(1).maxCoeff();
}

void Hash::addEdge(const Eigen::Vector2d& vi, const Eigen::Vector2d& vj,
    const Eigen::Vector2d& ui, const Eigen::Vector2d& uj, const int index)
{
    static Eigen::Vector2d edge_min, edge_max;
    calculate_edge_extents(vi, vj, ui, uj, edge_min, edge_max);
    this->addElement(edge_min, edge_max, index + 1); // Edges have a positive id
}

void Hash::addEdges(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i edges)
{
    for (int i = 0; i < edges.rows(); i++) {
        this->addEdge(vertices.row(edges(i, 0)), vertices.row(edges(i, 1)),
            displacements.row(edges(i, 1)), displacements.row(edges(i, 1)), i);
    }
}

void Hash::addElement(Eigen::Vector2d xmin, Eigen::Vector2d xmax, int id)
{
    AABBi aabbi = makeAABBi(xmin, xmax);
    this->add(aabbi, id);
}

void Hash::getVertexEdgePairs(const Eigen::MatrixX2i& edges, Candidates& hits)
{
    hits.clear();

    std::vector<int> edge_ids;
    edge_ids.reserve(edges.rows() / 10);
    std::vector<int> vertex_ids;
    vertex_ids.reserve(edges.rows() / 5);

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
            for (const int& vertex_id : vertex_ids) {
                for (const int& edge_id : edge_ids) {
                    if (edges(edge_id, 0) != vertex_id
                        && edges(edge_id, 1) != vertex_id) {
                        hits.insert(std::make_pair(vertex_id, edge_id));
                    }
                }
            }

            edge_ids.clear();
            vertex_ids.clear();
        }
    }
}

AABBi Hash::makeAABBi(Eigen::Vector2d mn, Eigen::Vector2d mx)
{
    AABBi aabbi;
    aabbi.m_min = Point2i(int(mn.x() / m_cellSize), int(mn.y() / m_cellSize));
    aabbi.m_max = Point2i(int(std::ceil(mx.x() / m_cellSize)),
        int(std::ceil(mx.y() / m_cellSize)));
    return aabbi;
}

void Hash::sort() { std::sort(m_hash.begin(), m_hash.end()); }

void Hash::add(AABBi aabbi, int id)
{
    for (int x = aabbi.getMin().x(); x <= aabbi.getMax().x(); ++x) {
        for (int y = aabbi.getMin().y(); y <= aabbi.getMax().y(); ++y) {
            add(Point2i(x, y), id);
        }
    }
}

void Hash::add(Point2i p, int id) { m_hash.push_back(HashItem(hash(p), id)); }

int Hash::hash(Point2i p) { return p.y() * m_gridSize + p.x(); }

void Hash::clear() { m_hash.clear(); }

HashItem& Hash::get(unsigned int i) { return m_hash[i]; }

} // namespace ccd
