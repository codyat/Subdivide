// CS130 Fall 2013: Computer Graphics
// point3d.h
//
// This file does not need to be modified
/////////////////////////////////////////
#ifndef __POINT2D_H__
#define __POINT2D_H__

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

//the index of each vertex in vertices
static int vertex_index = 0;

//the 3D point and location of the vertex
struct Point3D
{
	double x;
	double y;
	double z;
	
	Point3D() : x(0.0), y(0.0), z(0.0) {}
	Point3D(const double & nx, const double & ny, const double & nz) : x(nx), y(ny), z(nz) {}
	Point3D operator+(const Point3D & rhs) const { return Point3D(x + rhs.x, y + rhs.y, z + rhs.z); }
	Point3D operator-(const Point3D & rhs) const { return Point3D(x - rhs.x, y - rhs.y, z - rhs.z); }
	Point3D operator*(double val) const { return Point3D(x * val, y * val, z * val); }
	Point3D operator/(double val) const { return Point3D(x / val, y / val, z / val); }
	Point3D& operator+=(const Point3D & rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Point3D& operator-=(const Point3D & rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Point3D& operator*=(double val) { x *= val; y *= val; z *= val; return *this; }
	Point3D& operator/=(double val) { x /= val; y /= val; z /= val; return *this; }
	Point3D& operator=(const Point3D &rhs) {x = rhs.x; y = rhs.y; z = rhs.z;}
	
	float magnitude() const
	{ 
		return sqrt(x * x + y * y + z * z); 
	}
	void normalize()
	{ 
		*this /= magnitude(); 
	}
	float dot(const Point3D & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}
	Point3D cross(const Point3D & rhs) const
	{
		return Point3D(y * rhs.z - z * rhs.y,
					z * rhs.x - x * rhs.z,
					x * rhs.y - y * rhs.x);
	}
};

//a vertex consists of a 3D location, an index, and a vector of adjacent vertices
struct Vertex
{
	Point3D p;
	int index;
	vector<int> adj;
	
	Vertex() {}
	
	Vertex(const Vertex &rhs)
	{
		p.x = rhs.p.x;
		p.y = rhs.p.y;
		p.z = rhs.p.z;
		
		index = rhs.index;
		
		for(int i = 0; i < rhs.adj.size(); i++)	adj.push_back(rhs.adj[i]);
	}
	
	bool operator==(const Vertex &rhs)
	{
		return ((p.x == rhs.p.x) && (p.y == rhs.p.y) && (p.z == rhs.p.z));
	}
	
	Vertex& operator=(const Vertex &rhs)
	{
		if(this->p.x == rhs.p.x && this->p.y == rhs.p.y && this->p.z == rhs.p.z) 
			return *this; 
		
		this->p.x = rhs.p.x;
		this->p.y = rhs.p.y;
		this->p.z = rhs.p.z;
		
		this->index = rhs.index;
		this->adj.clear();
		
		for(int i = 0; i < rhs.adj.size(); i++)	this->adj.push_back(rhs.adj[i]);
		
		return *this;
	}
};

//the vector of vertices
vector<Vertex> vertices;

//an edge consist of p1, p2 (endpoints of the edge) p3, p4 (neighboring vertices) and a midpoint
//an edge can complete itself given p1 and p2
//
//unions calculates which vertices are common between the adjacency list of p1 and p2 and stores 
//those as p3 and p4
//
//midpoint calculates the 3d location for the midpoint of the edge as well as stating that the
//endpoints will immediately be adjacent to the midpoint (obviously)
struct Edge
{
	int p1;
	int p2;
	int p3;
	int p4;
	Vertex m;
	
	Edge() : p1(0), p2(0), p3(0), p4(0) {}
	Edge(int x, int y) : p1(x), p2(y), p3(0), p4(0) {}
	
	void midpoint()
	{
		m.p = vertices[p1].p*(.375) + vertices[p2].p*(.375) + vertices[p3].p*(.125) + vertices[p4].p*(.125);
		m.index = vertex_index++;
		m.adj.push_back(p1);
		m.adj.push_back(p2);
	}
	
	void unions()
	{
		int index = 0;
		int ind[2];
		for(int i = 0; i < vertices[p1].adj.size(); i++)
			for(int j = 0; j < vertices[p2].adj.size(); j++)
			{
				if(vertices[p1].adj[i] == vertices[p2].adj[j])
				{
					ind[index] = vertices[p1].adj[i];
					index++;
				}
			}
		p3 = ind[0];
		p4 = ind[1];
	}
	
	bool operator==(const Edge &rhs)
	{
		return (((p1 == rhs.p1) && (p2 == rhs.p2)) || ((p1 == rhs.p2) && (p2 == rhs.p1))); 
	}
};

//a face consists of 3 ints that index into the vector of edges 
struct Face
{
	int x;
	int y;
	int z;
	
	Face() : x(0), y(0), z(0) {}
	Face(const int & nx, const int & ny, const int & nz) : x(nx), y(ny), z(nz) {}
};

struct Ray
{
	Point3D R0;
	Point3D Rd;
	
	Ray()
	{
	  R0 = Point3D(0.0, 0.0, 0.0);
	  Rd = Point3D(0.0, 0.0, 1.0);
	}
};

#endif
