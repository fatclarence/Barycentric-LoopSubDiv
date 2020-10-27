#pragma once

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000

#include <math.h>
#include <utility>
#include <set>
#include <map>


using namespace std;

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version) { return (t << 3) | version; };
inline tIdx idx(OrTri ot) { return ot >> 3; };
inline int ver(OrTri ot) { return ot & 0b111; }; 
inline OrTri enext(OrTri ot) {
	int v = ver(ot);  return makeOrTri(idx(ot),
	                           v < 3 ? (v + 1) % 3 : 3 + ((v - 1) % 3)) ; };
inline OrTri sym(OrTri ot) { int v = ver(ot); return v < 3 ? ot + 3 : ot - 3; };

struct Vector3 {
	union {
		struct {
			double x, y, z;
		};
		double v[3];
	};

	Vector3() : x(0), y(0), z(0) {};

	Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

	double Magnitude() { return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); };
	void Normalize() {
		double mag = Magnitude();
		x /= mag;
		y /= mag;
		z /= mag;
	}
};

// Use of machine epsilon to compare floating-point values for equality
template<class T>
typename enable_if<!numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
	// the machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	return abs(x - y) <= numeric_limits<T>::epsilon() * abs(x + y) * ulp
		// unless the result is subnormal
		|| abs(x - y) < numeric_limits<T>::min();
}

class myObjType {
	int vcount = 0;
	int tcount = 0;

	double vlist[MAXV][3];   // vertices list
	int tlist[MAXT][3];      // triangle list
	int fnlist[MAXT][3];     // fnext list for future (not this assignment)
	double vnlist[MAXT][3];	 // vertex normal list
	double nlist[MAXT][3];   // storing triangle normals
	map<int, set<int>> adjVerticesOfVertex; // storing of neighbouring vertices of a vertex
	map<int, set<int>> adjFacesOfVertex; // storing of all triangles that share a vertex
	map<set<int>, set<int>> adjVerticesOfEdge; // store the left (and right if any) vertices of an edge

	// Storage variables for subdivisions
	int temp_tcount = 0;
	int temp_tlist[MAXT][3];

	double lmax[3];          // the maximum coordinates of x,y,z
	double lmin[3];          // the minimum coordinates of x,y,z

	int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
	int statMaxAngle[18]; 


public:
	myObjType() { vcount = 0; tcount = 0; };

	// Vector operations
	void normalize(double* v);
	double* crossProduct(double* v1, double* v2);
	double dotProduct(double* v1, double* v2);
	double* computeNormal(double* v1, double* v2, double* v3);
	double computeAngleBetweenVectors(double* v1, double* v2);

	Vector3 sumOfVectors(Vector3 v1, Vector3 v2);
	Vector3 scaleVector(double scale, Vector3 v);
	Vector3 doublesToVector3(double* v);
	double* vector3ToDoubles(Vector3 vec);

	void readFile(char* filename);  // assumming file contains a manifold
	void writeFile(char* filename);  
	void draw(bool isSmoothShade, bool isHighlightEdges);
	void highlightEdge();
    void computeStat();
	int org(OrTri ot);
	int dest(OrTri ot);
	void instantiateAdjacencyAttr();
	int computeComponents();
	void computeVertexNormals();
	void computeTriangleNormals();
	pair<int, int> getVertexPair(OrTri ot);
	bool checkOrientation(OrTri ot1, OrTri ot2);
	bool orientTriangles();

	// Subdivision Helper Functions
	Vector3 getInteriorOddVertex(int v1, int v2, int vleft, int vright);
	Vector3 getBoundaryOddVertex(int v1, int v2);
	
	pair<bool, int> addVertexToVlist(Vector3 vertex, double v_list[MAXV][3], int curr_vcount);
	void addTriangleToTlist(int temp_tlist[MAXT][3], int temp_tcount, int* vertices);
	void loopSubdivision();
	void barycentricSubdivision();
};

