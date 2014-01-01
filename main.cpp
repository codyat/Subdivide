// Name: Cody Troyer
// Quarter, Year: Fall, 2013
// Project 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include "point3d.h"
#include <fstream>

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const double PI = 3.14159265358979323846;
const double PI_FRACT = (PI)/90;

vector<Face> faces;
vector<Edge> edges;
ifstream in;
Point3D center;
Point3D light(800.0, 800.0, -6000.0);
Point3D camera(400.0, 400.0, 400);

int bmax_x, bmin_x, bmax_y, bmin_y;

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(double x, double y, float r = 0.0, float g = 0.5, float b = 0.5)
{
	glBegin(GL_POINTS);
	glColor3f(r,g,b);
	glVertex2f(x,y);
	glEnd();
}
//renders a line from one point to another regardless of slope
void DDA(double x1, double y1, double x2, double y2)
{

	double x = x2 - x1 + 0.0;
	double y = y2 - y1 + 0.0;
	double m = y / x;
	
	if((x2 < x1 && y2 < y1) || (x2 < x1 && y2 > y1))
	{
		double temp = x2;
		x2 = x1;
		x1 = temp;
		
		temp = y2;
		y2 = y1;
		y1 = temp;
	}
	if(m == 0 && x2 < x1)
	{
		double temp = x2;
		x2 = x1;
		x1 = temp;
		
		temp = y2;
		y2 = y1;
		y1 = temp;
	}
	if(m <= 0 && m > -1) // -1 > m > 0
	{
		double y = y1 + 0.0;
		for(double x = x1; x <= x2; ++x)
		{
			y = y + m;
			renderPixel(x, y);
		}
	}
	else if(m >= 1) // m > 1
	{
		double x = x1 + 0.0;
		for(double y = y1; y <= y2; ++y)
		{
			x = x + 1/m;
			renderPixel(x, y);
		}
	}
	else if(m > 0 && m < 1) // 0 < m < 1
	{
		double y = y1 + 0.0;
		for(double x = x1; x <= x2; ++x)
		{
			y = y + m;
			renderPixel(x, y);
		}
	}
	else // m < -1
	{
		double x = x2 + 0.0;
		for(double y = y2; y <= y1; ++y)
		{
			x = x + 1/m;
			renderPixel(x, y);
		}
	}
}
//determines if an edge is already in the vector of edges
//  if it is, it returns the index of the edge in the vector
//  if not, return the size of the vector
int is_edge_in(Edge a, const vector<Edge>&v)
{
	int i = 0;
	if(v.size() == 0) return i;
	else
		for(i = 0; i < v.size(); i++)
			if(a == v[i]) 
				return i;
	return i;
}
//literally orientates the triangle correctly
//the way the code os written, provides 3 edges to a face. and 2 endpoints 
//to each edge so 1 of 4 things can happen:
//  both edge 1 and 3 can be oriented correctly
//    we do nothing for this case
//  edge 1 is oriented incorrectly
//    swap endpoints of edge 1
//  edge 3 is oriented incorrectly
//    swap endpoints of edge 3
//  both edge 1 and 3 are oriented incorrectly
//    swap enpoints of both edge 1 and 3
//
//at this point we havent looked at edge 2 at all but we know 1 and 3 are 
//correct, so we can compare edge 2 to either edge 1 or 3 and have:
//  edge 2 is oriented correctly
//    do nothing for this case
//  edge 2 is oriented incorrectly
//    swap the endpoints of edge 2
void eliminate_ambiguity(const Face &f, Edge &t1, Edge &t2, Edge &t3)
{
    //both edges are correct, so keep the same data for the temp edges 
	if(edges[f.x].p1 == edges[f.z].p2)
	{
		t1 = edges[f.x];
		t3 = edges[f.z];
	}
    //edge 1 is incorrect, so swap its endpoints
	else if(edges[f.x].p2 == edges[f.z].p2)
	{
		t1.p1 = edges[f.x].p2;
		t1.p2 = edges[f.x].p1;
		
		t3 = edges[f.z];
	}
    //edge 3 is incorrect, so swap its endpoints
	else if(edges[f.x].p1 == edges[f.z].p1)
	{
		t1 = edges[f.x];
		
		t3.p1 = edges[f.z].p2;
		t3.p2 = edges[f.z].p1;
	}
    //both are incorect, so swap their endpoints
	else
	{
		t1.p1 = edges[f.x].p2;
		t1.p2 = edges[f.x].p1;
		
		t3.p1 = edges[f.z].p2;
		t3.p2 = edges[f.z].p1;
	}
	//edge 2 is correct, so leave it alone
	if(t3.p1 == edges[f.y].p2)
	{
		t2 = edges[f.y];
	}
    //edge 2 is incorrect, so swap its endpoints
	else 
	{
		t2.p1 = edges[f.y].p2;
		t2.p2 = edges[f.y].p1;
	}
}
//takes in a face along with a vector of edges and a vector of faces
//the purpose of this function is to split a face into 9 sub edges and 4 sub faces
void split_face(Face f, vector<Edge>& ev, vector<Face>& fv)
{
	Face face;
	Edge e, temp1, temp2, temp3;
	int tmp;
	
	eliminate_ambiguity(f, temp1, temp2, temp3);
	
	//make edges that consist of subface 1 and add it to their respective vectors
    
    //make the edge out of the given data for subface 1
	e.p1 = temp1.p1;
	e.p2 = edges[f.x].m.index;
    //check if the edge already exists in our mesh
	tmp = is_edge_in(e, ev);
    //if not, add it to the vector
	if(tmp == ev.size()) ev.push_back(e);
    //this edge is the first edge of the first face
	face.x = tmp;
	
    //repeat for the next edge of subface 1
	e.p1 = edges[f.x].m.index;
	e.p2 = edges[f.z].m.index;
    //check if the edge already exists in our mesh
	tmp = is_edge_in(e, ev);
	//if not, add it to the vector
    if(tmp == ev.size()) ev.push_back(e);
    //this edge is the second edge of the first face
    face.y = tmp;
	
    //repeat for the next edge of subface 1
	e.p1 = edges[f.z].m.index;
	e.p2 = temp3.p2;
    //check if the edge already exists in our mesh
	tmp = is_edge_in(e, ev);
    //if not, add it to the vector
    if(tmp == ev.size()) ev.push_back(e);
    //this edge is the third edge of the first face
    face.z = tmp;
	//sub face is complete, so add it to the vector of faces
	fv.push_back(face);
    
    ////////////////////////////////////////////////////////////
    //      repeat this process for the other 3 subfaces      //
	////////////////////////////////////////////////////////////
    
    //make edges that consist of face 2 and add it to their respective vectors
	e.p1 = temp2.p1;
	e.p2 = edges[f.y].m.index;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.x = tmp;

	e.p1 = edges[f.x].m.index;
	e.p2 = edges[f.y].m.index;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.y = tmp;

	e.p1 = edges[f.x].m.index;
	e.p2 = temp1.p2;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.z = tmp;

	fv.push_back(face);

	//make edges that consist of face 3 and add it to their respective vectors
	e.p1 = temp3.p1;
	e.p2 = edges[f.z].m.index;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.x = tmp;

	e.p1 = edges[f.y].m.index;
	e.p2 = edges[f.z].m.index;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.y = tmp;

	e.p1 = edges[f.y].m.index;
	e.p2 = temp2.p2;
	tmp = is_edge_in(e, ev);
	if(tmp == ev.size()) ev.push_back(e);
	face.z = tmp;

	fv.push_back(face);

	//make edges that consist of face 4 and add it to their respective vectors
	e.p1 = edges[f.x].m.index;
	e.p2 = edges[f.y].m.index;
	face.x = is_edge_in(e, ev);

	e.p1 = edges[f.y].m.index;
	e.p2 = edges[f.z].m.index;
	face.y = is_edge_in(e, ev);

	e.p1 = edges[f.z].m.index;
	e.p2 = edges[f.x].m.index;
	face.z = is_edge_in(e, ev);

	fv.push_back(face);
}
//applies weights on the "old" vertices
//the old_size represents the number of vertices in the previous mesh before 
//subdividing
void apply_weights(int old_size)
{
    //number of vertices adjacent to the current vertex
	int n;
    //how much the adjacent vertices influences the current vertex 
	double w;
    //for all the old vertices
	for(int i = 0; i < old_size; i++)
	{
		n = vertices[i].adj.size();
        //determine the proper weight
		if(n == 3) w = (3.0/16.0);
		else w = (3.0/48.0);
		
		Point3D newp;
        //traverse the adjacent vertices with the weights and apply it to 
        //the current vertex
		newp += vertices[i].p * (1 - n*w);
		for(int j = 0; j < vertices[i].adj.size(); j++)
			newp += vertices[vertices[i].adj[j]].p * w;
		//update the x, y, and z coordinates of the current vertex
		vertices[i].p.x = newp.x;
		vertices[i].p.y = newp.y;
		vertices[i].p.z = newp.z;
	}
}

//this function will do a single subdivide of the mesh
void subdivide()
{
  //creates containers to store the new vertices, edges, and faces
	//copies old vertices to the new vertices
  vector<Vertex> newv = vertices;
	vector<Edge> newe;
	vector<Face> newf;
	//////////////////////////////////////////////////////////////////////////////
	//                  midpoint section                                        //
	//////////////////////////////////////////////////////////////////////////////
	int old_size = vertices.size();
    //clears the adjacency list for all vertices in new vertices
	for(int i = 0; i < newv.size(); i++) newv[i].adj.clear();
    //adds the correct adjecent vertices based on the new subdivided mesh 
	for(int i = 0; i < edges.size(); i++)
	{
    //the midpoint of an edge will be adjacent to the endpoints of that edge 
    //and vice versa 
		newv[edges[i].p1].adj.push_back(edges[i].m.index);
		newv[edges[i].p2].adj.push_back(edges[i].m.index);
	}
    //adds the midpoints to the new vector of vertices
	for(int i = 0; i < edges.size(); i++) newv.push_back(edges[i].m);
	//a midpoint will also be adjacent to the other midpoints it shares a face with
    //so add them
    for(int i = 0; i < faces.size(); i++)
	{
        //midpoint on edge 1 is adjacent to midpoint on edge 2 and 3
		newv[edges[faces[i].x].m.index].adj.push_back(edges[faces[i].y].m.index);
		newv[edges[faces[i].x].m.index].adj.push_back(edges[faces[i].z].m.index);
		
        //midpoint on edge 2 is adjacent to midpoints on edge 1 and 3
		newv[edges[faces[i].y].m.index].adj.push_back(edges[faces[i].x].m.index);
		newv[edges[faces[i].y].m.index].adj.push_back(edges[faces[i].z].m.index);
		
        //midpoint on edge 3 is adjacent to midpoints on edge 1 and 2
		newv[edges[faces[i].z].m.index].adj.push_back(edges[faces[i].x].m.index);
		newv[edges[faces[i].z].m.index].adj.push_back(edges[faces[i].y].m.index);
	}
	//replace the old vertices with the new one
	vertices.swap(newv);
	//                             vertices are now updated
	
    //subdivide each individual face of the current mesh
	for(int i = 0; i < faces.size(); i++)
		split_face(faces[i], newe, newf);
	//replace the old edges with the new one
	edges.swap(newe);
    //replace the old faces with the new one
	faces.swap(newf);
	//                             edges and faces are now updated
	
	//////////////////////////////////////////////////////////////////////////////
	//                  weight section                                          //
	//////////////////////////////////////////////////////////////////////////////	
	//apply the weights of the old vertices
    apply_weights(old_size);
	//////////////////////////////////////////////////////////////////////////////
}

//this function will compute the vertices that are in common between an edges 
//endpoints adjacent vertices. then perform the midpoint of those edges. then 
//subdivide the mesh
void complete_edges()
{
    //for all edges
	for(int i = 0; i < edges.size(); i++)
	{ 
        //calculate vertices 3 and 4 for that edge
		edges[i].unions();	 
        //calculate the midpoint of that edge
		edges[i].midpoint();
	}
}

Ray ray;

float r,s,t;

Point3D eb;
Point3D ec;

int p1, p2, p3;

//this function calculates the Beta of the barycentric coordinate
float beta(int tri)
{
	p1 = edges[faces[tri].x].p1;
	p2 = edges[faces[tri].x].p2;
	if(edges[faces[tri].y].p1 == p1 || edges[faces[tri].y].p1 == p2) 
		p3 = edges[faces[tri].y].p2;
	else
		p3 = edges[faces[tri].y].p1;
	
	
	eb = vertices[p2].p - vertices[p1].p;
	ec = vertices[p3].p - vertices[p1].p;
	Point3D a = vertices[p1].p;
	
	float detA = (((eb.x * ec.y * ray.Rd.z) + 
								 (eb.z * ec.x * ray.Rd.y) +
								 (eb.y * ec.z * ray.Rd.x)) -
								((eb.z * ec.y * ray.Rd.x) + 
								 (eb.x * ec.z * ray.Rd.y) +
								 (eb.y * ec.x * ray.Rd.z)));
		
	float detB = ((((a.x - ray.R0.x) * -1 * ec.y * ray.Rd.z) + 
								 ((a.y - ray.R0.y) * -1 * ec.z * ray.Rd.x) +
								 ((a.z - ray.R0.z) * -1 * ec.x * ray.Rd.y)) -
								(((a.z - ray.R0.z) * -1 * ec.y * ray.Rd.x) + 
								 ((a.y - ray.R0.y) * -1 * ec.x * ray.Rd.z) +
								 ((a.x - ray.R0.x) * -1 * ec.z * ray.Rd.y)));
	
	return (detB/detA);
}

//this function calculates the Gamma of the barycentric coordinate
float gamma(int tri)
{
	p1 = edges[faces[tri].x].p1;
	p2 = edges[faces[tri].x].p2;
	if(edges[faces[tri].y].p1 == p1 || edges[faces[tri].y].p1 == p2) 
		p3 = edges[faces[tri].y].p2;
	else
		p3 = edges[faces[tri].y].p1;
	
	
	eb = vertices[p2].p - vertices[p1].p;
	ec = vertices[p3].p - vertices[p1].p;
	Point3D a = vertices[p1].p;
	
	float detA = (((eb.x * ec.y * ray.Rd.z) + 
								 (eb.z * ec.x * ray.Rd.y) +
								 (eb.y * ec.z * ray.Rd.x)) -
								((eb.z * ec.y * ray.Rd.x) + 
								 (eb.x * ec.z * ray.Rd.y) +
								 (eb.y * ec.x * ray.Rd.z)));
		
	float detC = (((-1 * eb.x * (a.y - ray.R0.y) * ray.Rd.z) + 
								 (-1 * eb.y * (a.z - ray.R0.z) * ray.Rd.x) +
								 (-1 * eb.z * (a.x - ray.R0.x) * ray.Rd.y)) -
								((-1 * eb.z * (a.y - ray.R0.y) * ray.Rd.x) + 
								 (-1 * eb.y * (a.x - ray.R0.x) * ray.Rd.z) +
								 (-1 * eb.x * (a.z - ray.R0.z) * ray.Rd.y)));
	
	return (detC/detA);
}

//this function calculates the t of the barycentric coordinate
float funt(int tri)
{
	p1 = edges[faces[tri].x].p1;
	p2 = edges[faces[tri].x].p2;
	if(edges[faces[tri].y].p1 == p1 || edges[faces[tri].y].p1 == p2) 
		p3 = edges[faces[tri].y].p2;
	else
		p3 = edges[faces[tri].y].p1;
	
	
	eb = vertices[p2].p - vertices[p1].p;
	ec = vertices[p3].p - vertices[p1].p;
	Point3D a = vertices[p1].p;
	
	float detA = (((eb.x * ec.y * ray.Rd.z) + 
								 (eb.z * ec.x * ray.Rd.y) +
								 (eb.y * ec.z * ray.Rd.x)) -
								((eb.z * ec.y * ray.Rd.x) + 
								 (eb.x * ec.z * ray.Rd.y) +
								 (eb.y * ec.x * ray.Rd.z)));
	
	float detT = (((eb.x * ec.y * (a.z - ray.R0.z)) + 
								 (eb.y * ec.z * (a.x - ray.R0.x)) +
								 (eb.z * ec.x * (a.y - ray.R0.y))) -
								((eb.z * ec.y * (a.x - ray.R0.x)) + 
								 (eb.y * ec.x * (a.z - ray.R0.z)) +
								 (eb.x * ec.z * (a.y - ray.R0.y))));
	
	return (detT/detA);
}

//this function determines whether or not the ray is traced into the triangle
//  if it is, it returns the depth of the triangle
//  else, it returns -600000
float is_in_triangle(int tri)
{
	float valt = funt(tri);
	float b = beta(tri);
	float g = gamma(tri);
	
	if(valt >= 0)
		if(b >= 0 && b <= 1 && g >= 0 && g <= 1)
			if(b + g <= 1)
				return valt;
	return -600000.0;
}
/*
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// R(t) = Ro + Rdt, t[0, infinity)                                            //
//                                                                            //
// B(r,s) = A + (B-A)r + (C-A)s                                               //
//                                                                            //
// Rdt - (B-A)r - (C-A)s = A-Ro                                               //
// Rdt -  (eb)r -  (ec)s = A-Ro                                               //
//                                                                            //
// r is BETA, s is GAMMA                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

*/

//this function determines which face is infront of which face
//the second part colors the object based on the amont of light it recieves
void intersection()
{	
	//shadow part
	/*
	 * bool shadows[800][800];
	ray.R0 = light;
	for(int m = bmin_x; m < bmax_x; m++)
	{
		for(int n = bmin_y; n < bmax_y; n++)
		{
			for(int i = 0; i < faces.size(); i++)
			{
				if(is_in_triangle(i) != -600000) 
					shadows[m][n] = true;
				else 
					shadows[m][n] = false;
			}
		}
	}
	cerr << "done" << endl;
	for(int m = 0; m < WINDOW_WIDTH; m++)
		for(int n = 0; n < WINDOW_HEIGHT; n++)
			if(!shadows[m][n])
				renderPixel(m, n, 0.0, 0.3, 0.0);
	cerr << "done2" << endl;
	ray.R0.z = 0.0; 
	ray.Rd = Point3D(0.0, 0.0, 1.0);
	*/
	for(int m = bmin_x; m <= bmax_x; m++)
	{ 
		ray.R0.x = m;
		for(int n = bmin_y; n <= bmax_y; n++)
		{
			ray.R0.y = n;
			float min = 9999999999.0;
			int index = -1;
			for(int i = 0; i < faces.size(); i++)
			{
				float temp = is_in_triangle(i);
				if(min > temp && temp != -600000.0)
				{
					min = temp;
					index = i;
				}	
			}
			if(index != -1)
			{
				p1 = edges[faces[index].x].p1;
				p2 = edges[faces[index].x].p2;
				if(edges[faces[index].y].p1 == p1 || edges[faces[index].y].p1 == p2) 
					p3 = edges[faces[index].y].p2;
				else
					p3 = edges[faces[index].y].p1;
				
				////////////////////////////////////////////////////////////////
				//
				// L = kd * Ld * (I [dot] n) + ks * Ls * (v [dot] r)^a + ka * La
				//
				////////////////////////////////////////////////////////////////
				
				//the color valuse for the diffuse, specular, and ambient lights 
				Point3D diffuse(.8, .4, .2);
				Point3D specular(.8, .4, .2);
				Point3D ambient(.8, .4, .2);
				
				//the shininess coefficient
				double a = 20.0;
				
				//2 edges of a face, used to calculate the normal of the face
				Point3D U(vertices[p3].p.x - vertices[p1].p.x, vertices[p3].p.y - vertices[p1].p.y, vertices[p3].p.z - vertices[p1].p.z);
				Point3D W(vertices[p2].p.x - vertices[p1].p.x, vertices[p2].p.y - vertices[p1].p.y, vertices[p2].p.z - vertices[p1].p.z);
				
				//creates the ilumination, normal, and viewer vectors of the current face
				Point3D I(light.x + m, light.y + n, light.z + min);
				Point3D N = U.cross(W);
				Point3D V(0,0,1);
				
				//creates the halfway vector
				Point3D h = I + V;
				//normalizes all the vectors being used
				h.normalize();
				I.normalize();
				N.normalize();
				
				//flips all negative normals
				if(N.z < 0) N = N * -1;
				//flips all normals, due to bug in the code......
				N *= -1;
				
				double temp = h.dot(N);
				
				//updates the diffuse and specular lighting
				diffuse *= I.dot(N);
				specular *= pow(temp, a);
								
				//sums up all three lights at the pixel
				Point3D phong = ambient + diffuse + specular;
				
				//colors the current pixel with the according color
				renderPixel(m, n, phong.x, phong.y, phong.z);
			}
		}
	}
}

void render_image()
{
	//////////////////////////////////////////////////////////////////////
	// 									pseudo code
	//////////////////////////////////////////////////////////////////////
	//
	//  for each pixel j
	//    determine closest intersection with objects in the scene	
	// 		if no hits
	// 			color pixel j background color
	// 		else
	//			illumination j for =0;
	// 			for each light k
	// 				compute shadow ray intersection
	// 				if hit
	//					illumination j += ambient k
	//				else
	//					find illumination k and add to illumination j
	// 			end light
	//		end hit
	//	end pixel
	// 
	//////////////////////////////////////////////////////////////////////
	
	//calculates the binding box of the object
	bmax_x = -99999;
	bmax_y = -99999;
	bmin_x = 99999;
	bmin_y = 99999;
	
	for(int i = 0; i < vertices.size(); i++)
	{ 
		if(vertices[i].p.x < double(bmin_x)) bmin_x = int(vertices[i].p.x);
		if(vertices[i].p.x > double(bmax_x)) bmax_x = int(vertices[i].p.x);
		if(vertices[i].p.y < double(bmin_y)) bmin_y = int(vertices[i].p.y);
		if(vertices[i].p.y > double(bmax_y)) bmax_y = int(vertices[i].p.y);
  }
  //makes sure the boundaries of the binding box are within the view pane
  if(bmin_x < 0) bmin_x = 0;
  if(bmax_x > 800) bmax_x = 800;
  if(bmin_y < 0) bmin_y = 0;
  if(bmax_y > 800) bmax_y = 800;
  
  //continues to the ray tracing
	intersection();
}

int mode = 0;
double tx, ty, tz;
//if the mouse is clicked, subdivide
void glKeyboard(unsigned char key, int x, int y)
{	
	if(key == 'r') //rotate
	{
  	mode = 1; 
  }
  else if(key == 't') //translate
	{
	  mode = 2;
  }
  else if(key == 'l') //subdivision level
	{
    //subdivide the mesh
		subdivide();	
    mode = 0;
  }
  else if(key == 'i') //render image
	{
		mode = -1;
  }
  if(mode > 0)
  {
		if(key == 'a') //left, move everything right
		{
	    if(mode == 1)
	    {
				//rotates each vertex and alters its x, y, and z values according 
				//to the rotation about which axis
				for(int i = 0; i < vertices.size(); i++) 
				{
					vertices[i].p -= center;
					tx = vertices[i].p.x;
					tz = vertices[i].p.z;
					vertices[i].p.x = tx * cos(PI_FRACT) + tz * sin(PI_FRACT);
					vertices[i].p.z = -1 * tx * sin(PI_FRACT) + tz * cos(PI_FRACT);
					vertices[i].p += center;
				}
			}
			else 
			{
				//increases the x values of all the vertices and centroid
				for(int i = 0; i < vertices.size(); i++) vertices[i].p.x+=5;
				center.x+=5;
			}
		}
		else if(key == 's') //down, move everything up
  	{
			if(mode == 1)
	    {
				//rotates each vertex and alters its x, y, and z values according 
				//to the rotation about which axis
				for(int i = 0; i < vertices.size(); i++) 
				{
					vertices[i].p -= center;
					ty = vertices[i].p.y;
					tz = vertices[i].p.z;
					vertices[i].p.y = ty * cos(-PI_FRACT) - tz * sin(-PI_FRACT);
					vertices[i].p.z = ty * sin(-PI_FRACT) + tz * cos(-PI_FRACT);
					vertices[i].p += center;
				}
			}
			else 
			{
				//increases the x values of all the vertices and centroid
				for(int i = 0; i < vertices.size(); i++) vertices[i].p.y+=5;
				center.y+=5;
			}
    }
    else if(key == 'd') //right, move everything left
	  {
			if(mode == 1)
	    {
				//rotates each vertex and alters its x, y, and z values according 
				//to the rotation about which axis
				for(int i = 0; i < vertices.size(); i++) 
				{
					vertices[i].p -= center;
					tx = vertices[i].p.x;
					tz = vertices[i].p.z;
					vertices[i].p.x = tx * cos(-PI_FRACT) + tz * sin(-PI_FRACT);
					vertices[i].p.z = -1 * tx * sin(-PI_FRACT) + tz * cos(-PI_FRACT);
					vertices[i].p += center;
				}
			}
			else 
			{
				//increases the x values of all the vertices and centroid
				for(int i = 0; i < vertices.size(); i++) vertices[i].p.x-=5;
				center.x-=5;
			}
    }
    else if(key == 'w') //up, move everything down
  	{
			if(mode == 1)
	    {
				//rotates each vertex and alters its x, y, and z values according 
				//to the rotation about which axis
				for(int i = 0; i < vertices.size(); i++) 
				{
					vertices[i].p -= center;
					ty = vertices[i].p.y;
					tz = vertices[i].p.z;
					vertices[i].p.y = ty * cos(PI_FRACT) - tz * sin(PI_FRACT);
					vertices[i].p.z = ty * sin(PI_FRACT) + tz * cos(PI_FRACT);
					vertices[i].p += center;
				}
			}
			else 
			{
				//increases the x values of all the vertices and centroid
				for(int i = 0; i < vertices.size(); i++) vertices[i].p.y-=5;
				center.y-=5;
			}
    }
	}
	glutPostRedisplay();
}

//Output function to OpenGL Buffer
void GL_render()
{
  glClear(GL_COLOR_BUFFER_BIT);
  
  //renders the image and not the mesh
  if(mode == -1)
	  render_image();
	else //draw each and every edge of the mesh
		for(int i = 0; i < edges.size(); i++)
      DDA(vertices[edges[i].p1].p.x, 
          vertices[edges[i].p1].p.y, 
          vertices[edges[i].p2].p.x, 
          vertices[edges[i].p2].p.y);
  glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    
	glutCreateWindow("CS 130 - ctroy001 Project 1");

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// For the purposes of this lab, this is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
	glutDisplayFunc(GL_render);
}

//determines if the index of the vertex is already in the vector of vertices
//if it is, returns true
//if not, returns false
bool is_int_in(int a, const vector<int>&v)
{
	if(v.size() == 0) return false;
	else
		for(int i = 0; i < v.size(); i++)
			if(a == v[i]) 
				return true;
	return false;
}

//reads in the original mesh from a file
void read_in()
{
    //number of vertices and number of faces in the mesh
	int v, f;
	in >> v >> f;
    
    //x, y, and z coordinates
	int x, y, z;
    //temporary point
	Point3D pt;
    //temporary vertex
	Vertex vert;
	
	 //read in all the x, y, and z coordinates and create a vertex
	for(int i = 0; i < v; i++)
	{
		in >> x >> y >> z;
		
		pt.x = double(8*x);
		pt.y = double(8*y);
		pt.z = double(8*z);
		vert.p = pt;
		vert.index = vertex_index++;
        //push the vertex into the vector
		vertices.push_back(vert);
	}
  
  //calculates the centeroid
  for(int i = 0; i < vertices.size(); i++)
  {
		center.x += vertices[i].p.x;
		center.y += vertices[i].p.y;
		center.z += vertices[i].p.z;
	}
	
	center.x /= v;
	center.y /= v;
	center.z /= v;
	
    //temporary vector that holds faces
	vector<Face> temp;
    //completes the adjacency list for all the vertices
	for(int i = 0; i < f; i++)
	{
		in >> x >> y >> z;
		
		temp.push_back(Face(x,y,z));

		if(!is_int_in(y, vertices[x].adj)) 
			vertices[x].adj.push_back(y);
		if(!is_int_in(z, vertices[x].adj)) 
			vertices[x].adj.push_back(z);
		
		if(!is_int_in(x, vertices[y].adj)) 
			vertices[y].adj.push_back(x);
		if(!is_int_in(z, vertices[y].adj)) 
			vertices[y].adj.push_back(z);
		
		if(!is_int_in(x, vertices[z].adj)) 
			vertices[z].adj.push_back(x);
		if(!is_int_in(y, vertices[z].adj)) 
			vertices[z].adj.push_back(y);
	}
    
    //a temporary edge
	Edge e;
    //uses the vector of faces above to create all the edges and faces of the mesh
	for(int k = 0; k < temp.size(); k++)
	{
		Face f;
		int tmp;
		
        //create and edge with the data given for the face
		e.p1 = temp[k].x;
		e.p2 = temp[k].y;
        //check if the edge is already in the edge vector
		tmp = is_edge_in(e, edges);
        //if not, add it
		if(tmp == edges.size()) edges.push_back(e);
		//give the edge as the first edge of the face
        f.x = tmp;
		
		//create and edge with the data given for the face
        e.p1 = temp[k].y;
		e.p2 = temp[k].z;
		//check if the edge is already in the edge vector
		tmp = is_edge_in(e, edges);
		//if not, add it
		if(tmp == edges.size()) edges.push_back(e);
		//give the edge as the second edge of the face
        f.y = tmp;
		
		//create and edge with the data given for the face
        e.p1 = temp[k].z;
		e.p2 = temp[k].x;
		//check if the edge is already in the edge vector
		tmp = is_edge_in(e, edges);
		//if not, add it
		if(tmp == edges.size()) edges.push_back(e);
		//give the edge as the third edge of the face
        f.z = tmp;
		
        //pushes the completed face into the vector of faces
		faces.push_back(f);
	}
	//input from file is complete
	in.close();
	complete_edges();
	
}

int main(int argc, char** argv)
{	
	in.open(argv[1]);
	GLInit(&argc, argv);
	glutKeyboardFunc(glKeyboard);
	read_in();
	glutDisplayFunc(GL_render);
	glutMainLoop();

	return 0;
}
