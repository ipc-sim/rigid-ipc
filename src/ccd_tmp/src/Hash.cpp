// Hash.cpp
//

#include "Hash.h"
#include <iostream>

using namespace std;
// using namespace __gnu_cxx;
using namespace Eigen;

void Hash::resize(Vector3d mn, Vector3d mx, double cellSize)
{
    clear();
    m_cellSize = cellSize * 2.0;
    m_domainMin = mn;
    m_domainMax = mx;
    m_gridSize = std::max((mx[0] - mn[0]) / m_cellSize,
        std::max((mx[1] - mn[1]) / m_cellSize,
            (mx[2] - mn[2]) / m_cellSize));
}

void Hash::addElement(Vector3d xmin, Vector3d xmax, int id)
{
    AABBi aabbi = makeAABBi(xmin, xmax);
    add(aabbi, id);
}

void Hash::getVertexEdgePairs(unsigned int* e, Candidates& hits)
{
    hits.clear();

    // NOTE: If larger meshes are used this should be increased, or
    // the memory should be allocated dynamically
    //
    static size_t edge[10000];
    static size_t vs[10000];

    size_t eiSize = 0;
    size_t viSize = 0;

    // Sorted all they (key,value) pairs, where key is the hash key,
    // and value is the element index
    //
    sort();

    unsigned size = m_hash.size();
    assert(size != 0);

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    //
    for (unsigned i = 0; i < size; ++i) {
        int& currId = get(i).id;
        int& currH = get(i).key;

        // read this element
        if (currId > 0) // Edge elements have positive id
        {
            edge[eiSize++] = currId - 1;
            assert(eiSize < 10000);
        } else // Vertex elements have negative id
        {
            vs[viSize++] = -currId - 1;
            assert(viSize < 10000);
        }

        // CLOSE BUCKET: if this is the last element, or next is from another bucket
        if ((i == size - 1) || currH != get(i + 1).key)
        {
            // We are closing the bucket (key entry), so tally up
            // all vertex-edge pairs encountered in the bucket that just ended
            for (size_t vi = 0; vi < viSize; ++vi)
                for (size_t ei = 0; ei < eiSize; ++ei) {
                    // std::cout << "vs " << vs[vi] << " edge " << edge[ei] << std::endl;
                    if (e[2 * edge[ei]] != vs[vi] && e[2 * edge[ei] + 1] != vs[vi]) {
                        hits.push_back(make_pair(vs[vi], edge[ei]));
                    }
                }

            eiSize = 0;
            viSize = 0;
        }
    }
}

AABBi Hash::makeAABBi(Vector3d mn, Vector3d mx)
{
    AABBi aabbi;
    aabbi.m_min = Point3i(mn[0] / m_cellSize, mn[1] / m_cellSize, mn[2] / m_cellSize);
    aabbi.m_max = Point3i(mx[0] / m_cellSize, mx[1] / m_cellSize, mx[2] / m_cellSize);
    return aabbi;
}

void Hash::sort()
{
    std::sort(m_hash.begin(), m_hash.end());
}

void Hash::add(AABBi aabbi, int id)
{
    //std::cout<<"inside hash: id in hash obj: "<<id<<std::endl;
    for (int x = aabbi.getMin().x(); x <= aabbi.getMax().x(); ++x)
        for (int y = aabbi.getMin().y(); y <= aabbi.getMax().y(); ++y)
            for (int z = aabbi.getMin().z(); z <= aabbi.getMax().z(); ++z) {
                add(Point3i(x, y, z), id);
            }
}

void Hash::add(Point3i p, int id)
{
    m_hash.push_back(HashItem(hash(p), id));
}

int Hash::hash(Point3i p)
{
    return p.z() * m_gridSize * m_gridSize + p.y() * m_gridSize + p.x();
}

void Hash::clear()
{
    m_hash.clear();
}

HashItem& Hash::get(unsigned int i)
{
    return m_hash[i];
}
