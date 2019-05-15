// Curve.cpp
//

#include "Curve.h"
#include "InterferenceVolume.h"

#include <vector>
#include <iostream>
#include <ext/hash_map>

#include <Eigen/Core>
#include "rpoly.h"
#include "Hash.h"
#include <unordered_set>
using namespace std;
using namespace __gnu_cxx;
using namespace Eigen;



namespace IAGM
{

void Curve::findAndRemoveInterference()
{
	size_t nbrIters = 0;

	Intersections is;
	while (getCurveIntersections(is) && nbrIters < m_maxNbrIters)
	{
		InterferenceVolumes ivs;
		modelInterference(is, ivs);
		removeInterference(ivs);

		copyControlConfiguration();
		updateSlaveMesh();

		nbrIters++;
		is.clear();
	}

	if (nbrIters == m_maxNbrIters)
	{
		// TODO: Practically speaking, it might be reasonable to cap
		// the number of iterations to something reasonable. If this
		// threshold is hit, then rollback to the start configuration,
		// which is guaranteed to be safe (since the start configuration
		// contains to intersections, by hypothesis).
		//
		cerr << "Warning: maximum number of iterations exceeded! "
		     << "Mesh may contain intersections." << endl;
	}
}

void Curve::modelInterferenceGrouped(Intersections &is, InterferenceVolumes &ivs){
    // @orig: At this point one could partition the intersections into disjoint
	// regions, in order to create one interference volume per region.
	// DONE, by @libin
	vector<Intersections> groups;
	groupIntersections(is,groups);
	int mycount = 0;
    
	for (vector<Intersections>::iterator isItr=groups.begin(); isItr!=groups.end(); ++isItr)
	{
        modelInterference(*isItr, ivs);
    }
}
void Curve::modelInterference(Intersections &is, InterferenceVolumes &ivs)
{
	// We tag each vertex in the mesh with the earliest intersection in its
    // barycentric region that it is involved in.
    //
    
    size_t idx = 0;
    hash_map<size_t, size_t> h;
//    std::cout << "=================================================\n";
    for (IntersectionsIterator iItr=is.begin();iItr!=is.end(); ++iItr, ++idx)
    {
        if (h.find(iItr->getVertex()) == h.end() ||
                iItr->getTime() < is[h[iItr->getVertex()]].getTime())
            h[iItr->getVertex()] = idx;

        if (iItr->getAlpha() <= 0.5 &&
            (h.find(iItr->getEdgeVertex1()) == h.end() ||
                iItr->getTime() < is[h[iItr->getEdgeVertex1()]].getTime()))
            h[iItr->getEdgeVertex1()] = idx;
        
        if (iItr->getAlpha() >= 0.5 &&
            (h.find(iItr->getEdgeVertex2()) == h.end() ||
                iItr->getTime() < is[h[iItr->getEdgeVertex2()]].getTime()))
            h[iItr->getEdgeVertex2()] = idx;

//        std::cout << "node="<< iItr->getVertex()
//                  << " edge_0=" << iItr->getEdgeVertex1()
//                  << " edge_1=" << iItr->getEdgeVertex2()
//                  << "\n";
    }
//    for(auto it = h.begin(); it != h.end(); ++it)
//    {
//        std::cout << "vtx="<< it->first << " idx=" << it->second << "\n";
//    }

    // Create interference volume
    ivs.push_back(InterferenceVolume(this, is, h));
    
}


void Curve::groupIntersections(Intersections &inputInters, vector<Intersections> &intersections)
{
	intersections.clear();

	// First put every intersection in its own list
	//
	std::list<Intersections> intersectionsLists;
	std::list<unordered_set<unsigned int> > vertexLists;
	for (IntersectionsIterator iItr=inputInters.begin(); iItr!=inputInters.end(); ++iItr)
	{

		Intersections intersc;
		intersc.push_back(*iItr);
		intersectionsLists.push_back(intersc);

		unordered_set<unsigned int> s;
		s.insert(iItr->getVertex());
		s.insert(iItr->getEdgeVertex1());
		s.insert(iItr->getEdgeVertex2());


		// To be done 
		// insert neighboring vertex which are connected to one of the three vertices by some edge
		int verts[3];
		verts[0] = iItr->getVertex();
		verts[1] = iItr->getEdgeVertex1();
		verts[2] = iItr->getEdgeVertex2();
		
		vertexLists.push_back(s);
	}
	

	// Merge intersections with common vertices
	//
	while (!intersectionsLists.empty())
	{
		// Start with first intersection "zone"
		//
		Intersections intersc = intersectionsLists.front();
		intersectionsLists.pop_front();
		unordered_set<unsigned int> si = vertexLists.front();
		vertexLists.pop_front();

		int s = intersectionsLists.size();
		bool hits = false;

		// Loop through other zones
		//
		for (int i=0; i<s; ++i)
		{

			Intersections n = intersectionsLists.front();
			intersectionsLists.pop_front();
			unordered_set<unsigned int> next = vertexLists.front();
			vertexLists.pop_front();

			// Search for common elements
			//
			bool common = false;
			for (unordered_set<unsigned int>::iterator sItr=next.begin(); sItr!=next.end(); ++sItr)
			{

				if (si.find(*sItr) != si.end())
				{
					// They share a vertex
					//
					common = true;
					break;
				}

			}

			// If there are elements in common, then take the union of the two sets
			//
			if(common)
			{
				// Merge
				//
				si.insert(next.begin(),next.end());
				intersc.insert(intersc.end(), n.begin(), n.end());

				// We need to keep track of new zones for later, so signify that there's a new one
				//
				hits = true;
			}
			else
			{
				// Otherwise put it back on the queue
				//
				vertexLists.push_back(next);
				intersectionsLists.push_back(n);
			}
		}

		// If the search returned no hits, then put this set on the final list
		//
		if (!hits)
			intersections.push_back(intersc);
		else
		{
			intersectionsLists.push_back(intersc);
			vertexLists.push_back(si);
		}
	}

}


void Curve::removeInterference(InterferenceVolumes &ivs)
{
	double *ctrl=0; // control mesh vertices
	getControlConfiguration(ctrl);

	set<size_t> movable;
	getMovableControlVertices(movable);

	for (InterferenceVolumesIterator ivItr=ivs.begin(); ivItr!=ivs.end(); ++ivItr)
	{
		ivItr->computeHandleGradient(this, movable);

		double p = ivItr->getPressure();
		for (set<size_t>::iterator sItr=movable.begin(); sItr!=movable.end(); ++sItr)
		{
			for (size_t i=0; i<2; ++i)
				ctrl[2*(*sItr)+i] += p * ivItr->getHandleGradients()[2*(*sItr)+i];
		}
	}
}

bool Curve::getCurveIntersections(Intersections &is)
{
	double *x0=0; // Start configuration
	double *x1=0; // End configuration
    size_t *xGroup = 0; // vertices' group
	size_t nbrVerts = getNbrVertices();
	getStartConfiguration(x0);
	getEndConfiguration(x1);
    getGroup(xGroup);

	unsigned int *e=0; // Edge i consists of vertices e[2*i] and e[2*i+1]
	size_t nbrEdges;
	getEdgeIndices(nbrEdges, e);

	static Hash grid;

	// Resize grid
	//
	Vector3d mn = Vector3d( 1e100,  1e100, 0.0);
	Vector3d mx = Vector3d(-1e100, -1e100, 1.0);
	for (size_t i=0; i<nbrVerts; ++i)
	{
		mn[0] = min(x0[2*i+0], min(x1[2*i+0], mn[0]));
		mn[1] = min(x0[2*i+1], min(x1[2*i+1], mn[1]));
		mx[0] = max(x0[2*i+0], max(x1[2*i+0], mx[0]));
		mx[1] = max(x0[2*i+1], max(x1[2*i+1], mx[1]));
	}

	// TODO: For increased performance cells should probably be sized based
	// on average edge length (or something)
	//
	double cellSize = getCellSize();
	grid.resize(mn, mx, cellSize);

	// Insert points into grid
	//
	double *xStart = &x0[0];
	double *xEnd   = &x1[0];
	for (int i=0; i<(int)nbrVerts; ++i)
	{
		Vector3d xMin, xMax;
		for (size_t j=0; j<2; ++j)
		{
			xMin[j] = min(xStart[j], xEnd[j]) - getMinimumSeparation() - 1.0e-6;
			xMax[j] = max(xStart[j], xEnd[j]) + getMinimumSeparation() + 1.0e-6;
			//xMin[j] = min(xStart[j], xEnd[j]) - getMinimumSeparation();
			//xMax[j] = max(xStart[j], xEnd[j]) + getMinimumSeparation();
		}
		xMin[2] = -getMinimumSeparation();
		xMax[2] =  getMinimumSeparation();

		// To differentiate between vertices and edges, vertex indices
		// will always be negative (and offset by 1, to prevent ambiguity
		// for the 0-th element)
		//
		grid.addElement(xMin, xMax, -i-1);

		xStart += 2;
		xEnd   += 2;
	}

	// Insert edges into grid
	//
	for (size_t i=0; i<nbrEdges; ++i)
	{
		double *x10 = &x0[2*e[2*i+0]];
		double *x11 = &x1[2*e[2*i+0]];
		double *x20 = &x0[2*e[2*i+1]];
		double *x21 = &x1[2*e[2*i+1]];
		
		Vector3d xMin, xMax;
		for (size_t j=0; j<2; ++j)
		{
			xMin[j] = min(x10[j], min(x11[j], min(x20[j], x21[j]))) - getMinimumSeparation() - 1.0e-6;
			xMax[j] = max(x10[j], max(x11[j], max(x20[j], x21[j]))) + getMinimumSeparation() + 1.0e-6;
		}
		xMin[2] = -getMinimumSeparation();
		xMax[2] =  getMinimumSeparation();

		grid.addElement(xMin, xMax, i+1);
	}

	// Query potential intersections
	//
	Candidates hits;
	grid.getVertexEdgePairs(e, hits);


	// Low-level intersection test
	//
	is.clear();
	int mycount_is = 0;
	for (CandidatesIterator cItr=hits.begin(); cItr!=hits.end(); ++cItr)
	{
		double s, t;
        //  @libin: avoid self intersections
        if(m_ignore_self_intersections && (xGroup[2*cItr->first] == xGroup[2*e[2*cItr->second+0]])){
            continue;
        }

		if (getIntersection(&x0[2*cItr->first],         &x1[2*cItr->first],
							&x0[2*e[2*cItr->second+0]], &x1[2*e[2*cItr->second+0]],
							&x0[2*e[2*cItr->second+1]], &x1[2*e[2*cItr->second+1]], s, t))
			{
				is.push_back(Intersection(cItr->first, e[2*cItr->second], e[2*cItr->second+1], s, t));
				mycount_is = mycount_is + 1;
			}
	}
	return !is.empty();
}

bool Curve::getIntersection(double *x00, double *x01,
							double *x10, double *x11, double *x20, double *x21, double &s, double &t)
{

	Vector2d x1(x00[0], x00[1]);
	Vector2d x2(x10[0], x10[1]);
	Vector2d x3(x20[0], x20[1]);
	
	Vector2d v1(x01[0]-x00[0], x01[1]-x00[1]);
	Vector2d v2(x11[0]-x10[0], x11[1]-x10[1]);
	Vector2d v3(x21[0]-x20[0], x21[1]-x20[1]);

	Vector2d x1x2 = x1 - x2;
	Vector2d v1v2 = v1 - v2;
	Vector2d x3x2 = x3 - x2;
	Vector2d v3v2 = v3 - v2;
	
	swap(x3x2[0], x3x2[1]); x3x2[1] *= -1.0;
	swap(v3v2[0], v3v2[1]); v3v2[1] *= -1.0;

	double a0 = x1x2.dot(x3x2);
	double a1 = x1x2.dot(v3v2) + x3x2.dot(v1v2);
	double a2 = v1v2.dot(v3v2);

	vector<double> coeffs;
  	if (getMinimumSeparation() == 0.0)
  	{
    		coeffs.resize(3);
    		coeffs[0] = a2;
    		coeffs[1] = a1;
    		coeffs[2] = a0;
  	}
  	else
  	{
        	double h = getMinimumSeparation();
            // DEBUG
            // std::cout << "h " << h << std::endl;
    		coeffs.resize(5);

    		coeffs[0] = a2 * a2;
    		coeffs[1] = 2.0 * a1 * a2;
	    	coeffs[2] = a1 * a1 + 2.0 * a0 * a2;
	    	coeffs[3] = 2.0 * a0 * a1;
	    	coeffs[4] = a0 * a0;
	
	    	coeffs[2] -= h * h * v3v2.dot(v3v2);
	    	coeffs[3] -= h * h * 2.0 * x3x2.dot(v3v2);
	    	coeffs[4] -= h * h * x3x2.dot(x3x2);
  	}

	// TODO: Coefficient pruning tests could speed this up
	//

	int degree = coeffs.size() - 1;
  	while (coeffs[coeffs.size() - 1 - degree] == 0.0 && degree > 0)
		degree--;

	int cnt = 0;
	double rl[4], im[4];
	RootFinder rf;
  	if (degree == 0)
	{
		cnt = 0;
	}
	else
	{
		// Find roots of polynomial
		//
		cnt = rf.rpoly(&coeffs[coeffs.size() - 1 - degree], degree, rl, im);	
	}

	// Try vertex-segment two half circles collision
	// TODO now it's checking the entire circle, should
	// check half circle, checking the entire circle
	// may cause problems.
	Vector2d x1x3 = x1 - x3;
	Vector2d v1v3 = v1 - v3;
	double minSep = getMinimumSeparation();

	// Test x1 collides wtih half circle of x3
	double coeffs_circle1[3];
	coeffs_circle1[0] = v1v3.dot(v1v3);
	coeffs_circle1[1] = 2*x1x3.dot(v1v3);
	coeffs_circle1[2] = x1x3.dot(x1x3) - minSep*minSep;
	degree = 2;
	
	if (coeffs_circle1[0] == 0.0)
		degree = 1;
	
	double rl_circle1[2], im_circle1[2];
	int cnt_circle1 = rf.rpoly(&coeffs_circle1[2-degree], degree, rl_circle1, im_circle1);

	// Test x1 collides with half circle of x2
	double coeffs_circle2[3];
	coeffs_circle2[0] = v1v2.dot(v1v2);
	coeffs_circle2[1] = 2*x1x2.dot(v1v2);
	coeffs_circle2[2] = x1x2.dot(x1x2) - minSep*minSep;
	degree = 2;
	
	if (coeffs_circle2[0] == 0.0)
		degree = 1;
	
	double rl_circle2[2], im_circle2[2];
	int cnt_circle2 = rf.rpoly(&coeffs_circle2[2-degree], degree, rl_circle2, im_circle2);

	//cnt_circle1 = 0;
	//cnt_circle2 = 0;
	int cnt_total = cnt + cnt_circle1 + cnt_circle2;
	double rl_total[8], im_total[8];
	int root_count = 0;
	for (int i = 0; i<cnt; ++i)
	{
		rl_total[root_count] = rl[i];
		im_total[root_count] = im[i];
		root_count = root_count + 1;
	}
	
	for (int i = 0; i<cnt_circle1; ++i)
	{
		rl_total[root_count] = rl_circle1[i];
		im_total[root_count] = im_circle1[i];
		root_count = root_count + 1;
	}

	for (int i = 0; i<cnt_circle2; ++i)
	{
		rl_total[root_count] = rl_circle2[i];
		im_total[root_count] = im_circle2[i];
		root_count = root_count + 1;
	}


	double min_t = 10000;
	double min_s = 10000;
	double min_d = 10000;
	
	Vector2d pe = x1 + v1;
	Vector2d ae = x2 + v2;
	Vector2d be = x3 + v3;

  	Vector2d abe = be - ae;
	Vector2d ape = pe - ae;

	double se = abe.dot(ape) / abe.dot(abe);
	if (se < 0.00)   se = 0.0;
	if (se > 1.00)   se = 1.0;
	
	double de = (pe - (ae + abe * se)).squaredNorm() - getMinimumSeparation() * getMinimumSeparation();

	Vector2d ps = x1;
	Vector2d as = x2;
	Vector2d bs = x3;

  	Vector2d abs = bs - as;
	Vector2d aps = ps - as;

	double ss = abs.dot(aps) / abs.dot(abs);
	if (ss < 0.00)   ss = 0.0;
	if (ss > 1.00)   ss = 1.0;
	
	double ds = (ps - (as + abs * ss)).squaredNorm() - getMinimumSeparation() * getMinimumSeparation();
	// //DEBUG
	// if(ds<-1e-10)
	//   std::cout << "ds " << ds <<std::endl;
	
	for (int i=0; i<cnt_total; ++i)
	{
		//if (im_total[i] == 0.0 && rl_total[i] >= 0.0 && rl_total[i] <= 1.0)
		if ( (std::abs(im_total[i]) < 1e-16) && (rl_total[i] >= -1e-16) && (rl_total[i] <= (1.0+1e-16)))
		{
			t = rl_total[i];
                        if(t<0) t = 0;
			if(t>1) t = 1;

			// Check rl[i] for intersection
			//
			Vector2d p = x1 + v1 * t;
			Vector2d a = x2 + v2 * t;
			Vector2d b = x3 + v3 * t;

			Vector2d ab = b - a;
			Vector2d ap = p - a;

			s = ab.dot(ap) / ab.dot(ab);

			//std::cout<<"s: "<<s<<std::endl;
			// Clamp to edge endpoints
			//
			if (s < 0.00)   s = 0.0;
			if (s > 1.00)   s = 1.0;


			// NOTE: If the code misses intersections, it's probably because of
			// this tolerance. There is a Eurographics 2012 paper which attempts to
			// address this, but for now it can be tricky to choose a reliable epsilon.
			//
			double d = (p - (a + ab * s)).squaredNorm() - getMinimumSeparation() * getMinimumSeparation();
            // DEBUG distance between nodes.
            // std::cout << "d " << d << std::endl;

			// safeguard test if t is close to 1, check the ending position
			if ((d < m_tol) || (de <= 0))
			{
				// safeguard test if t is close to 0, check the start position 

				// To ensure a conservative interference volume, slightly enlarge
				// it by returning an earlier time
				// t = max(0.0, t - 1e-6);
				if(t < min_t)
				{
					min_t = t;
					min_s = s;
					min_d = d;
				}
			}
		}
	}

	if(min_t >= 0.0 && min_t <= 1.0)
	{
		t = min_t;
		s = min_s;
		return true;
	
	}

	return false;
}


} // end namespace IAGM

