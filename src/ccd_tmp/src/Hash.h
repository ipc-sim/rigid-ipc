// Hash.h
//
// This is a slightly modified version of code written by
// Daniele Panozzo
//

#ifndef _HASH_H_
#define _HASH_H_

#include <vector>
#include <list>
#include <iomanip>

#include <Eigen/Core>

typedef std::pair<size_t, size_t> Candidate;
typedef std::vector<Candidate> Candidates;
typedef std::vector<Candidate>::iterator CandidatesIterator;



// An entry into the hash grid
//
class HashItem
{
public:
	HashItem(int k, int i)
	{
		key = k;
		id = i;
	}

	HashItem() {}

	int key;
	int id;

	bool operator<(const HashItem& hi) const
	{
		return key < hi.key;
	}
};

template<class T>
class AABBT
{
public:
	
	AABBT(const T& min, const T& max)
	{
		m_min = min;
		m_max = max;
	}
	
	AABBT()
	{
		m_min = T(0,0,0);
		m_max = T(0,0,0);
	}
	
	~AABBT() { }
	
	T& getMin()
	{ return m_min; }
	
	T& getMax()
	{ return m_max; }
	
	T m_min; 
	T m_max; 
};

typedef Eigen::Matrix<int,3,1> Point3i;
typedef Eigen::Matrix<double,3,1> Point3d;
typedef AABBT<Point3i> AABBi;
typedef AABBT<Point3d> AABBd;



class Hash
{
public:

	void resize(Eigen::Vector3d mn, Eigen::Vector3d mx, double cellSize);

	void addElement(Eigen::Vector3d xmin, Eigen::Vector3d xmax, int id);

	void getVertexEdgePairs(unsigned int *e, Candidates &hits);

protected:
	AABBi makeAABBi(Eigen::Vector3d mn, Eigen::Vector3d mx);

	void sort();

	void add(AABBi aabbi, int id);
	void add(Point3i p, int id);

	int hash(Point3i p);		
	void clear();
	
	HashItem& get(unsigned int i);

protected:
	double m_cellSize;
	int m_gridSize;
	Eigen::Vector3d m_domainMin;
	Eigen::Vector3d m_domainMax;

	std::vector<HashItem> m_hash;

};

#endif

