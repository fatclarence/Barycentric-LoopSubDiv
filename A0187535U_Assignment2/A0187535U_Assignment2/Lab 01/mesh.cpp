#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <iomanip>
#include <sstream>
#include <set>
#include <utility>
#include <algorithm>

using namespace std;

////////////////////////////////////////////////////////////
// VECTOR OPERATIONS
///////////////////////////////////////////////////////////


void myObjType::normalize(double* v)
{
	double length = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));

	for (int i = 0; i < 3; i++) {
		v[i] = v[i] / length;
	}
}

double* myObjType::crossProduct(double* v1, double* v2)
{
	double result[3] = {
		v1[1] * v2[2] - v1[2] * v2[1],
		-(v1[0] * v2[2] - v1[2] * v2[0]),
		v1[0] * v2[1] - v1[1] * v2[0]
	};

	return result;
}

double myObjType::dotProduct(double* v1, double* v2)
{
	normalize(v1);
	normalize(v2);
	double product = 0;

	for (int k = 0; k < 3; k++) {
		product = product + (v1[k] * v2[k]);
	}

	return product;
}

// Compute face normal
double* myObjType::computeNormal(double* v1, double* v2, double* v3)
{
	double v1_to_v2[3] = { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
	double v1_to_v3[3] = { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };

	normalize(v1_to_v2);
	normalize(v1_to_v3);

	double* result = crossProduct(v1_to_v2, v1_to_v3);
	return result;
}

double myObjType::computeAngleBetweenVectors(double* v1, double* v2)
{
	double dot_product_of_vectors = dotProduct(v1, v2);
	return acos(dot_product_of_vectors) * (180 / M_PI);
}

////////////////////////////////////////////////////////////
// END OF VECTOR OPERATIONS
///////////////////////////////////////////////////////////

void myObjType::draw(bool isSmoothShade, bool isHighlightEdges) {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);

	if (isHighlightEdges) {
		highlightEdge();
		glPopMatrix();
		return;
	}

	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);

		if (!isSmoothShade) {
			glNormal3dv(nlist[i]);
		}
		
		for (int j = 0; j < 3; j++) {
			if (isSmoothShade) {
				glNormal3dv(vnlist[tlist[i][j]]);
			}
			glVertex3dv(vlist[tlist[i][j]]);
		}
		glEnd();
	
	}
	glDisable(GL_LIGHTING);

	glPopMatrix();
}

void myObjType::highlightEdge() {
	glDisable(GL_LIGHTING);

	for (int i = 0; i <= tcount; i++) {
		for (int version = 0; version < 3; version++) {
			// Check if the triangle has any neighbours
			OrTri neighbour = fnlist[i][version];
			if (idx(neighbour) == 0) {
				// Get the edge vertex
				OrTri current = makeOrTri(i, version);
				pair<int, int> edge_vertex_pair = getVertexPair(current);
				glBegin(GL_LINES);
				glColor3f(1.0f, 0.0f, 0.0f);
				glVertex3dv(vlist[edge_vertex_pair.first]);
				glVertex3dv(vlist[edge_vertex_pair.second]);
				glEnd();
			}
		}
	}
}

void myObjType::writeFile(char* filename)
{
	// .OBJ file format
	if (strstr(filename, ".obj") != NULL) {

		ostringstream file_lines;

		// Write all vertices
		for (int i = 1; i <= vcount; i++) {
			file_lines << "v " 
					<< to_string(vlist[i][0]) << " " 
					<< to_string(vlist[i][1]) << " " 
					<< to_string(vlist[i][2]) << endl;
		}

		// Write all faces
		for (int i = 1; i <= tcount; i++) {
			file_lines << "f "
				<< to_string(tlist[i][0]) << " "
				<< to_string(tlist[i][1]) << " "
				<< to_string(tlist[i][2]) << endl;
		}
		ofstream outFile(filename);
		outFile << file_lines.str() << endl;

		outFile.close();
	}
	// .OFF file format
	else if (strstr(filename, ".off") != NULL) {

		ostringstream file_lines;

		file_lines << "OFF" << endl;

		file_lines << vcount << " " << tcount << endl;

		// Write all vertices
		for (int i = 1; i <= vcount; i++) {
			file_lines << to_string(vlist[i][0]) << " "
				<< to_string(vlist[i][1]) << " "
				<< to_string(vlist[i][2]) << endl;
		}

		// Write all faces
		// Note that in .off files, vertices start at index of 0
		int vertices_per_face = 3;
		for (int i = 1; i <= tcount; i++) {
			file_lines << vertices_per_face << " "
				<< to_string(tlist[i][0] - 1) << " "
				<< to_string(tlist[i][1] - 1) << " "
				<< to_string(tlist[i][2] - 1) << endl;
		}
		ofstream outFile(filename);
		outFile << file_lines.str();

		outFile.close();
	}
	else {
		cerr << "Error: File format " << filename << " not supported!";
		exit(1);
	}
}

void myObjType::readFile(char* filename)
{
	if (strstr(filename, ".obj") != NULL) {
		cout << "Opening " << filename << endl;
		ifstream inFile;
		inFile.open(filename);
		if (!inFile.is_open()) {
			cout << "We cannot find your file " << filename << endl;
			exit(1);
		}

		string line;
		int i, j;
		bool firstVertex = 1;
		double currCood;

		while (getline(inFile, line))
		{
			if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
			{
				if (line[0] == 'v')
				{
					vcount++;
					i = 1;
					const char* linec = line.data();
					for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
						while (linec[i] == ' ') i++;
						j = i;
						while (linec[j] != ' ') j++;
						currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str()); // Convert C String to a Double
						if (firstVertex)
							lmin[k] = lmax[k] = currCood;
						else {
							if (lmin[k] > currCood)
								lmin[k] = currCood;
							if (lmax[k] < currCood)
								lmax[k] = currCood;
						}
						i = j;
					}

					firstVertex = 0;
				}
				if (line[0] == 'f')
				{
					tcount++;
					i = 1;
					const char* linec = line.data();
					for (int k = 0; k < 3; k++) {
						while (linec[i] == ' ') i++;
						j = i;
						while (linec[j] != ' ' && linec[j] != '\\') j++;
						tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
						i = j;
						fnlist[tcount][k] = 0;
						while (linec[j] != ' ') j++;
					}

				}


			}
		}
	} else if (strstr(filename, ".off") != NULL) {
		cout << "Opening " << filename << endl;
		ifstream inFile;
		inFile.open(filename);
		if (!inFile.is_open()) {
			cout << "We cannot find your file " << filename << endl;
			exit(1);
		}

		string line;
		int i, j, number_of_vertices, number_of_faces, current_line = 0;
		bool firstVertex = 1;
		double currCood;

		while (getline(inFile, line))
		{
			// Check if first line is OFF
			if (current_line == 0) {
				if (line.substr(0, 3) != "OFF") {
					cerr << "Error: Corrupted .off file. Please check your input .off file again." << endl;
					exit(1);
				}
				current_line++;

			} else if (current_line == 1) {
				// Second line of OFF file contains the number of vertices and faces
				const char* linec = line.data();
				for (int j = 0; j < 2; j++) {
					i = 0;
					while (linec[i] != ' ') i++;
					
					// Save number of vertices
					number_of_vertices = atof(line.substr(0, i).c_str());

					// Save number of triangles
					number_of_faces = atof(line.substr(i + 1, line.length()).c_str());
				}
				current_line++;

			} else {
				// Storing vertices
				if (current_line <= number_of_vertices + 1) {
					vcount++;
					i = 0;
					const char* linec = line.data();
					for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
						while (linec[i] == ' ') i++;
						j = i;
						while (linec[j] != ' ') j++;
						// vlist starts from index 1
						currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str()); // Convert C String to a Double
						if (firstVertex)
							lmin[k] = lmax[k] = currCood;
						else {
							if (lmin[k] > currCood)
								lmin[k] = currCood;
							if (lmax[k] < currCood)
								lmax[k] = currCood;
						}
						i = j;
					}

					firstVertex = 0;
					current_line++;
					
				} 

				// Storing faces
				else if (current_line > number_of_vertices + 1 && current_line <= number_of_vertices + number_of_faces + 2) {
					tcount++;
					i = 1;
					const char* linec = line.data();
					for (int k = 0; k < 3; k++) {
						while (linec[i] == ' ') i++;
						j = i;
						while (linec[j] != ' ' && linec[j] != '\\') j++;
						// tlist starts at index 1
						tlist[tcount][k] = atof(line.substr(i, j - i).c_str()) + 1;
						i = j;
						fnlist[tcount][k] = 0;
						while (linec[j] != ' ') j++;
					}
					current_line++;
				}
			}
			
		}
	}

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;

    computeStat();
	instantiateAdjacencyAttr();
	computeTriangleNormals();
	computeVertexNormals();

	computeComponents();

	orientTriangles();

	cout << "Done" << endl;
}

void myObjType::computeStat()
{
	cout << "Computing statistics..." << endl;
	int i;
    double minAngle = 360;
    double maxAngle = 0;

	for (int j = 1; j <= tcount; j++) {
		int v1, v2, v3;
		v1 = tlist[j][0];
		v2 = tlist[j][1];
		v3 = tlist[j][2];

		double v1_to_v2[3] = { vlist[v2][0] - vlist[v1][0], vlist[v2][1] - vlist[v1][1], vlist[v2][2] - vlist[v1][2] };
		double v1_to_v3[3] = { vlist[v3][0] - vlist[v1][0], vlist[v3][1] - vlist[v1][1], vlist[v3][2] - vlist[v1][2] };
		double v2_to_v3[3] = { vlist[v3][0] - vlist[v2][0], vlist[v3][1] - vlist[v2][1], vlist[v3][2] - vlist[v2][2] };
		double v2_to_v1[3] = { -(vlist[v2][0] - vlist[v1][0]), -(vlist[v2][1] - vlist[v1][1]), -(vlist[v2][2] - vlist[v1][2]) };

		double angle1 = computeAngleBetweenVectors(v1_to_v3, v1_to_v2);
		double angle2 = computeAngleBetweenVectors(v2_to_v1, v2_to_v3);
		double angle3 = 180.0 - angle1 - angle2;

		double current_min = min(min(angle1, angle2), angle3);
		double current_max = max(max(angle1, angle2), angle3);

		minAngle = min(current_min, minAngle);
		maxAngle = max(current_max, maxAngle);

		statMinAngle[int(floor(current_min / 10))] += 1;
		statMaxAngle[int(floor(current_max / 10))] += 1;
	}

    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;

	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;
}

int myObjType::org(OrTri ot) 
{
	// Find triangle index of OrTri
	tIdx triangle_index = idx(ot);

	// Find version of this triangle
	int version = ver(ot);

	// Map of versions and their org vertices
	map<int, int> oriented_triangles_map = {
		{0, 0}, {1, 1}, {2, 2},
		{3, 1}, {4, 2}, {5, 0}
	};

	int origin = tlist[triangle_index][oriented_triangles_map[version]];

	return origin;
}

int myObjType::dest(OrTri ot)
{
	return org(sym(ot));
}

void myObjType::instantiateAdjacencyAttr()
{
	cout << "Instantiating Adjacency Attributes..." << endl;
	map<set<int>, set<OrTri>> edge_to_triangle;

	// Record the triangles that share one edge
	// Also record the triangles that surround a vertex
	for (int i = 1; i <= tcount; i++)
	{
		for (int version = 0; version < 3; version++)
		{
			int v0 = tlist[i][version];
			int v1 = tlist[i][(version + 1) % 3];
			int v2 = tlist[i][(version + 2) % 3];
			OrTri ot = makeOrTri(i, version);
			set<int> edge = { v0, v1 };

			edge_to_triangle[edge].insert(ot);
			
			// Tag an edge with the adjacent vertices (if any)
			if (v2 != 0) {
				adjVerticesOfEdge[edge].insert(v2);
			}

			// Tag this triangle's index to the current vertex
			adjFacesOfVertex[v0].insert(i);

			// Tag the other 2 vertex that form this triangle as neighbours of current vertex
			adjVerticesOfVertex[v0].insert(v1);
			adjVerticesOfVertex[v0].insert(v2);
		}
	}

	cout << "All neighbours for all vertices discovered" << endl;
	cout << "All adjacent vertices of all edges discovered" << endl;

	// Read back the map of triangles and edges to form fnlist
	for (int i = 1; i <= tcount; i++)
	{
		for (int version = 0; version < 3; version++)
		{
			int v0 = tlist[i][version];
			int v1 = tlist[i][(version + 1) % 3];

			set<int> edge = { v0, v1 };

			// 2 triangles recorded per edge
			set<OrTri> faces_sharing_an_edge = edge_to_triangle[edge];

			OrTri f0, f1;

			f0 = *next(faces_sharing_an_edge.begin(), 0);

			// If there is only 1 face that has this edge, the edge is at the boundary
			if (faces_sharing_an_edge.size() == 1)
				f1 = 0;
			else
				f1 = *next(faces_sharing_an_edge.begin(), 1);

			int fnext = i == idx(f0) ? f1 : f0;
			fnlist[i][version] = fnext;
		}
	}

	cout << "Fnlist instantiated." << endl;
}

void myObjType::computeVertexNormals() {

	cout << "Computing vertex normals" << endl;
	// Read back map of triangles to instantiate vnlist
	for (int i = 1; i <= vcount; i++) {

		// Get all triangles that surround current vertex
		set<int> surrounding_triangles = adjFacesOfVertex[i];

		Vector3 normal;

		for (int k = 0; k < surrounding_triangles.size(); k++) {
			int triangle_index = *next(surrounding_triangles.begin(), k);
			Vector3 v = doublesToVector3(nlist[triangle_index]);
			normal = sumOfVectors(normal, v);
		}

		normal.Normalize();

		vnlist[i][0] = normal.x;
		vnlist[i][1] = normal.y;
		vnlist[i][2] = normal.z;
	}

	cout << "Vertex normals computed" << endl;
}

void myObjType::computeTriangleNormals() {
	cout << "Compute triangle normals" << endl;

	for (int i = 1; i <= tcount; i++) {
		int v1, v2, v3;
		v1 = tlist[i][0];
		v2 = tlist[i][1];
		v3 = tlist[i][2];

		double* normalVector = computeNormal(vlist[v1], vlist[v2], vlist[v3]);

		nlist[i][0] = normalVector[0];
		nlist[i][1] = normalVector[1];
		nlist[i][2] = normalVector[2];
	}

	cout << "Computed triangle normals" << endl;
}


int myObjType::computeComponents()
{
	cout << "Computing components..." << endl;
	// Perform BFS
	int number_of_components = 0;
	queue<int> indices_to_traverse;
	set<int> traversed_indices;

	while (traversed_indices.size() < tcount)
	{
		int index_of_undiscovered_triangle = 0;
		for (int i = 1; i <= tcount; i++) {
			if (traversed_indices.find(i) == traversed_indices.end())
			{
				index_of_undiscovered_triangle = i;
				break;
			} else {
				index_of_undiscovered_triangle = 0;
			}
		}

		if (index_of_undiscovered_triangle == 0)
		{
			break;
		}
		else 
		{
			indices_to_traverse.push(index_of_undiscovered_triangle);
		}

		while (indices_to_traverse.size() != 0)
		{
			int curr_discovered_triangle = indices_to_traverse.front();
			indices_to_traverse.pop();
			
			// Find all the neighbouring triangles and add their indices into the queue
			for (OrTri& neighbour : fnlist[curr_discovered_triangle]) {
				if (idx(neighbour) != 0 && traversed_indices.find(curr_discovered_triangle) == traversed_indices.end()) {
					indices_to_traverse.push(idx(neighbour));
				}
			}

			traversed_indices.insert(curr_discovered_triangle);
		}

		// After traversal of all connected triangles, 1 component is found
		number_of_components++;
	}

	cout << "Number of Components: " << number_of_components << endl;
	return number_of_components;
}

pair<int, int> myObjType::getVertexPair(OrTri ot) {

	int version = ver(ot);
	int index = idx(ot);
	int v0, v1;

	if (version == 0) {
		v0 = tlist[index][0];
		v1 = tlist[index][1];
	}
	else if (version == 1) {
		v0 = tlist[index][1];
		v1 = tlist[index][2];
	}
	else {
		v0 = tlist[index][2];
		v1 = tlist[index][0];
	}

	return make_pair(v0, v1);
}

bool myObjType::checkOrientation(OrTri ot1, OrTri ot2) {
	pair<int, int> t1VertexPair = getVertexPair(ot1);
	
	pair<int, int> t2VertexPair = getVertexPair(ot2);

	// If the vertex pair are the same, yet ot1 != ot2 then ot1 and ot2 have different orientations
	return t1VertexPair == t2VertexPair; 
}

bool myObjType::orientTriangles() {
	cout << "Execute orient triangles..." << endl;
	// Perform BFS
	int number_of_oriented = 0;
	queue<int> indices_to_traverse;
	set<int> traversed_indices;

	while (traversed_indices.size() < tcount)
	{
		int index_of_undiscovered_triangle = 0;
		for (int i = 1; i <= tcount; i++) {
			if (traversed_indices.find(i) == traversed_indices.end())
			{
				index_of_undiscovered_triangle = i;
				break;
			}
			else {
				index_of_undiscovered_triangle = 0;
			}
		}

		if (index_of_undiscovered_triangle == 0)
		{
			break;
		}
		else
		{
			indices_to_traverse.push(index_of_undiscovered_triangle);
		}

		while (indices_to_traverse.size() != 0)
		{
			int curr_discovered_triangle = indices_to_traverse.front();
			indices_to_traverse.pop();

			// Find all the neighbouring triangles and add their indices into the queue
			for (int version = 0; version < 3; version++) {

				OrTri neighbour = fnlist[curr_discovered_triangle][version];
				OrTri current = makeOrTri(curr_discovered_triangle, version);
				bool is_different = checkOrientation(current, neighbour);
				int neighbour_index = idx(neighbour);
				if (traversed_indices.find(curr_discovered_triangle) == traversed_indices.end()) {

					indices_to_traverse.push(neighbour_index);

					if (is_different) {

						// Orientate triangles by reversing vertices
						int prev_vertex = tlist[neighbour_index][1];
						tlist[neighbour_index][1] = tlist[neighbour_index][2];
						tlist[neighbour_index][2] = prev_vertex;

						// Update fnlist: Version 2 will have the neighbours of version 0 of triangle
						OrTri prev = fnlist[neighbour_index][0];
						fnlist[neighbour_index][0] = fnlist[neighbour_index][2];
						fnlist[neighbour_index][2] = prev;

						number_of_oriented++;
					}
				}
				else {

				}
			}
			traversed_indices.insert(curr_discovered_triangle);
		}

	}

	cout << "Number of Triangles Oriented: " << number_of_oriented << endl;

	

	// Compute new normals and update nlist
	for (int i = 1; i <= tcount; i++) {
		int v1, v2, v3;
		v1 = tlist[i][0];
		v2 = tlist[i][1];
		v3 = tlist[i][2];

		double* normalVector = computeNormal(vlist[v1], vlist[v2], vlist[v3]);

		nlist[i][0] = normalVector[0];
		nlist[i][1] = normalVector[1];
		nlist[i][2] = normalVector[2];
	}

	instantiateAdjacencyAttr();
	computeTriangleNormals();
	computeVertexNormals();

	return true;
}


//////////////////////////////////////////////////////////////////
// LOOP SUBDIVISION CODE
/////////////////////////////////////////////////////////////////

Vector3 myObjType::sumOfVectors(Vector3 v1, Vector3 v2) {
	Vector3 sum;

	sum.x = v1.x + v2.x;
	sum.y = v1.y + v2.y;
	sum.z = v1.z + v2.z;

	return sum;
}

Vector3 myObjType::scaleVector(double scale, Vector3 v) {
	Vector3 result;

	result.x = v.x * scale;
	result.y = v.y * scale;
	result.z = v.z * scale;

	return result;
}

Vector3 myObjType::doublesToVector3(double* v) {
	Vector3 result;

	result.x = v[0];
	result.y = v[1];
	result.z = v[2];

	return result;
}

double* myObjType::vector3ToDoubles(Vector3 vec) {
	double v[3];

	v[0] = vec.x;
	v[1] = vec.y;
	v[2] = vec.z;
	
	return v;
}

/**
* Functions for Odd Loop Vertex Computation
*/

Vector3 myObjType::getInteriorOddVertex(int v1, int v2, int vleft, int vright) {

	Vector3 edge_vertex1, edge_vertex2, left_vertex, right_vertex;

	// Instantiate vectors
	edge_vertex1 = doublesToVector3(vlist[v1]);
	edge_vertex2 = doublesToVector3(vlist[v2]);
	left_vertex = doublesToVector3(vlist[vleft]);
	right_vertex = doublesToVector3(vlist[vright]);

	// Scale with weights proposed by Warren for odd vertices
	// Store intermediate results
	Vector3 intermediate1 = scaleVector(3. / 8., sumOfVectors(edge_vertex1, edge_vertex2));
	Vector3 intermediate2 = scaleVector(1. / 8., sumOfVectors(left_vertex, right_vertex));

	Vector3 new_vertex = sumOfVectors(intermediate1, intermediate2);
	return new_vertex;
}

Vector3 myObjType::getBoundaryOddVertex(int v1, int v2) {

	Vector3 edge_vertex1, edge_vertex2;

	edge_vertex1 = doublesToVector3(vlist[v1]);
	edge_vertex2 = doublesToVector3(vlist[v2]);

	Vector3 result = scaleVector(1. / 2., sumOfVectors(edge_vertex1, edge_vertex2));
	return result;
}

// Adds new vertex if it is not already into vlist.
pair<bool, int> myObjType::addVertexToVlist(Vector3 vertex, double v_list[MAXV][3], int curr_vcount) {
	for (int i = 1; i <= curr_vcount; i++) {
		if (almost_equal(vertex.x, v_list[i][0], 2)
			&& almost_equal(vertex.y, v_list[i][1], 2)
			&& almost_equal(vertex.z, v_list[i][2], 2))

			return make_pair(false, i);
	}

	v_list[curr_vcount + 1][0] = vertex.x;
	v_list[curr_vcount + 1][1] = vertex.y;
	v_list[curr_vcount + 1][2] = vertex.z;
	return make_pair(true, curr_vcount + 1);
}

// Adds a new triangle into the tlist.
void myObjType::addTriangleToTlist(int temp_tlist[MAXT][3], int temp_tcount, int* vertices) {
	temp_tlist[temp_tcount + 1][0] = vertices[0];
	temp_tlist[temp_tcount + 1][1] = vertices[1];
	temp_tlist[temp_tcount + 1][2] = vertices[2];
}

void myObjType::loopSubdivision() {

	cout << "Executing Loop Subdivision..." << endl;
	for (int i = 1; i <= tcount; i++) {
		vector<Vector3> odd_vertices;
		int even_vertices[3];

		for (int version = 0; version < 3; version++) {
			int edge_vidx1 = tlist[i][version];
			int edge_vidx2 = tlist[i][(version + 1) % 3];

			set<int> edge = { edge_vidx1, edge_vidx2 };

			set<int> adjacent_vertices = adjVerticesOfEdge[edge];

			Vector3 odd_vertex, even_vertex;

			if (adjacent_vertices.size() == 1) {
				// Boundary edge case
				odd_vertex = getBoundaryOddVertex(edge_vidx1, edge_vidx2);
			}
			else {
				// Interior edge case
				int vleft = *next(adjacent_vertices.begin(), 0);
				int vright = *next(adjacent_vertices.begin(), 1);
				odd_vertex = getInteriorOddVertex(edge_vidx1, edge_vidx2, vleft, vright);
			}

			odd_vertices.push_back(odd_vertex);

			even_vertices[version] = edge_vidx1;
		}

		vector<int> added_vertices;
		
		// Add odd vertices
		for (int k = 0; k < 3; k++) {
			pair<bool, int> result = addVertexToVlist(odd_vertices[k], vlist, vcount);
			if (result.first) {
				vcount++;
			}
			added_vertices.push_back(result.second);
		}

		// Add even vertices
		added_vertices.insert(added_vertices.end(), begin(even_vertices), end(even_vertices));

		int t1[3] = { added_vertices[3], added_vertices[0], added_vertices[2] };
		int t2[3] = { added_vertices[0], added_vertices[4], added_vertices[1] };
		int t3[3] = { added_vertices[2], added_vertices[0], added_vertices[1] };
		int t4[3] = { added_vertices[2], added_vertices[1], added_vertices[5] };

		int* triangles[4] = { t1, t2, t3, t4 };

		for (int k = 0; k < 4; k++) {
			addTriangleToTlist(temp_tlist, temp_tcount++, triangles[k]);
		}
 	}

	//copy(&temp_vlist[0][0], &temp_vlist[0][0] + (temp_vcount + 1) * 3, &vlist[0][0]);
	copy(&temp_tlist[0][0], &temp_tlist[0][0] + (temp_tcount + 1) * 3, &tlist[0][0]);

	tcount = temp_tcount;

	cout << "Complete loop subdivision" << endl;
	cout << "No. of vertices: " << vcount << endl;
	cout << "No. of triangles: " << tcount << endl;

	instantiateAdjacencyAttr();
	computeTriangleNormals();
	computeVertexNormals();

	orientTriangles();

	cout << "Done" << endl;

	// Reset all temporary variables
	temp_tcount = 0;
}


// Form 6 triangles using a centroid and add all the triangles into the tlist.
void myObjType::barycentricSubdivision() {
	cout << "Executing Barycentric Subdivision" << endl;

	for (int i = 1; i <= tcount; i++) {
		vector<Vector3> vertices_in_triangle;
		vector<int> vidx_in_triangle;
		Vector3 sum;

		for (int version = 0; version < 3; version++) {
			int edge_idx1 = tlist[i][version];
			int edge_idx2 = tlist[i][(version + 1) % 3];

			vidx_in_triangle.push_back(edge_idx1);

			Vector3 edge_vertex1 = doublesToVector3(vlist[edge_idx1]);
			Vector3 edge_vertex2 = doublesToVector3(vlist[edge_idx2]);

			Vector3 mid_point = scaleVector(1./2., sumOfVectors(edge_vertex1, edge_vertex2));

			vertices_in_triangle.push_back(edge_vertex1); //0, 2, 4 are vertices that are already in vlist
			vertices_in_triangle.push_back(mid_point);

			sum = sumOfVectors(sum, edge_vertex1);
		}

		Vector3 centroid = scaleVector(1. / 3., sum);
		vertices_in_triangle.push_back(centroid);


		for (int k = 0; k < vertices_in_triangle.size() - 1; k++) {
			// add the vertex index of all mid points into the vidx_in_triangle arr
			if (k % 2 != 0) {
				Vector3 v = vertices_in_triangle[k];
				pair<bool, int> result = addVertexToVlist(v, vlist, vcount);
				if (result.first) {
					vcount++;
				}
				vidx_in_triangle.push_back(result.second);
			}
		}

		// Add centroid into vidx_in_triangle arr
		Vector3 v = vertices_in_triangle[6];
		pair<bool, int> result = addVertexToVlist(v, vlist, vcount);
		if (result.first) {
			vcount++;
		}
		vidx_in_triangle.push_back(result.second);

		// Form 6 new triangles with the 7 vertices
		int t1[3] = { vidx_in_triangle[6], vidx_in_triangle[0], vidx_in_triangle[3] };
		int t2[3] = { vidx_in_triangle[6], vidx_in_triangle[3], vidx_in_triangle[1] };
		int t3[3] = { vidx_in_triangle[6], vidx_in_triangle[1], vidx_in_triangle[4] };
		int t4[3] = { vidx_in_triangle[6], vidx_in_triangle[4], vidx_in_triangle[2] };
		int t5[3] = { vidx_in_triangle[6], vidx_in_triangle[2], vidx_in_triangle[5] };
		int t6[3] = { vidx_in_triangle[6], vidx_in_triangle[5], vidx_in_triangle[0] };

		int* triangles[6] = { t1, t2, t3, t4, t5, t6 };

		for (int k = 0; k < 6; k++) {
			addTriangleToTlist(temp_tlist, temp_tcount++, triangles[k]);
		}

	}

	cout << "All new triangles formed" << endl;

	cout << "Adding triangles into tlist..." << endl;

	copy(&temp_tlist[0][0], &temp_tlist[0][0] + (temp_tcount + 1) * 3, &tlist[0][0]);

	tcount = temp_tcount;

	cout << "Complete barycentric subdivision" << endl;
	cout << "No. of vertices: " << vcount << endl;
	cout << "No. of triangles: " << tcount << endl;

	instantiateAdjacencyAttr();
	computeTriangleNormals();
	computeVertexNormals();

	orientTriangles();

	temp_tcount = 0;

	cout << "Done" << endl;
}