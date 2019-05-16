// Hash.cpp
//

#include "hash.hpp"
#include <iostream>

using namespace std;
// using namespace __gnu_cxx;
using namespace Eigen;

void Hash::resize(Vector2d mn, Vector2d mx, double cellSize)
{
    clear();
    m_cellSize = cellSize * 2.0;
    m_domainMin = mn;
    m_domainMax = mx;
    m_gridSize = std::max(mx.x() - mn.x(), mx.y() - mn.y()) / m_cellSize;
}

void Hash::addElement(Vector2d xmin, Vector2d xmax, int id)
{
    AABBi aabbi = makeAABBi(xmin, xmax);
    add(aabbi, id);
}

void Hash::getVertexEdgePairs(const Eigen::MatrixX2i& edges, Candidates& hits)
{
    hits.clear();

    std::vector<long> edge_ids;
    std::vector<long> vertex_ids;

    // Sorted all they (key,value) pairs, where key is the hash key,
    // and value is the element index
    //
    sort();

    // unsigned size = m_hash.size();
    // assert(size != 0);

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    //
    for (unsigned i = 0; i < m_hash.size(); ++i) {
        const int& currId = get(i).id;
        const int& currH = get(i).key;

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
            // for (size_t vi = 0; vi < viSize; ++vi) {
            for (const long& vertex_id : vertex_ids) {
                // for (size_t ei = 0; ei < eiSize; ++ei) {
                for (const long& edge_id : edge_ids) {
                    // std::cout << "vs " << vs[vi] << " edge " << edge[ei] <<
                    // std::endl;
                    if (edges(edge_id, 0) != vertex_id
                        && edges(edge_id, 1) != vertex_id) {
                        hits.push_back(make_pair(vertex_id, edge_id));
                    }
                }
            }

            edge_ids.clear();
            vertex_ids.clear();
        }
    }
}

AABBi Hash::makeAABBi(Vector2d mn, Vector2d mx)
{
    AABBi aabbi;
    aabbi.m_min = Point2i(mn.x(), mn.y()) / m_cellSize;
    aabbi.m_max = Point2i(mx.x(), mx.y()) / m_cellSize;
    return aabbi;
}

void Hash::sort() { std::sort(m_hash.begin(), m_hash.end()); }

void Hash::add(AABBi aabbi, int id)
{
    // std::cout<<"inside hash: id in hash obj: "<<id<<std::endl;
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
