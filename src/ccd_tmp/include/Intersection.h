// Intersection.h
//
// This class represents an intersection, which is a point in space-time
// when two geometric elements touch. This class can represent an intersection
// between a vertex-edge pair, a vertex-triangle pair, and an edge-edge pair
//

#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include <vector>



namespace IAGM
{

class Intersection
{
public:
    Intersection(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3,
                 double t, bool vf)
	 : m_v0(v0), m_v1(v1), m_v2(v2), m_v3(v3), m_time(t), m_vf(vf)
    { }
    
	Intersection(unsigned int v0, unsigned int v1, unsigned int v2, double s, double t)
	 : m_v0(v0), m_v1(v1), m_v2(v2), m_alpha(s), m_beta(-1.0), m_gamma(-1.0), m_time(t)
    { }

	Intersection()
	 : m_vf(true)
    { }

	bool isVF() const
	{ return m_vf; }

	bool isEE() const
	{ return !m_vf; }
	
	double getTime() const
	{ return m_time; }

	void set(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3)
	{ m_v0 = v0; m_v1 = v1; m_v2 = v2; m_v3 = v3; }
	
	void set(unsigned int v0, unsigned int v1, unsigned int v2)
	{ m_v0 = v0; m_v1 = v1; m_v2 = v2; }

	// Edge-edge collision accessors
	//
	unsigned int getFirstEdgeVertex1() const
	{ return m_v0; }
	
	unsigned int getFirstEdgeVertex2() const
	{ return m_v1; }
	
	unsigned int getSecondEdgeVertex1() const
	{ return m_v2; }
	
	unsigned int getSecondEdgeVertex2() const
	{ return m_v3; }

	// Vertex-triangle accessors
	//
    unsigned int getVertex() const
    { return m_v0; }
    
    unsigned int getTriangleVertex1() const
    { return m_v1; }
    
    unsigned int getTriangleVertex2() const
    { return m_v2; }
    
    unsigned int getTriangleVertex3() const
    { return m_v3; }

	// Vertex-edge accessors
	//
	unsigned int getEdgeVertex1() const
	{ return m_v1; }

	unsigned int getEdgeVertex2() const
	{ return m_v2; }

	// Parametrized value accessors
	//
	double getAlpha() const
	{ return m_alpha; }

	double getBeta() const
	{ return m_beta; }

	double getGamma() const
	{ return m_gamma; }

protected:
    unsigned int m_v0;
    unsigned int m_v1;
    unsigned int m_v2;
    unsigned int m_v3;

	double m_alpha;
	double m_beta;
	double m_gamma;

    double m_time;
	
	bool m_vf;

private:

};

typedef std::vector<Intersection> Intersections;
typedef std::vector<Intersection>::iterator IntersectionsIterator;

} // end namespace IAGM

#endif

