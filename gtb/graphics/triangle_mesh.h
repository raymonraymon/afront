
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/



#ifndef __TRIANGLE_MESH_H
#define __TRIANGLE_MESH_H

#include <gtb/common.hpp>
#include <gtb/graphics/model.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/box3.hpp>



GTB_BEGIN_NAMESPACE

template <class T>
class tTriangleMeshVertex {
public:
    typedef T value_type;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;

	tTriangleMeshVertex() : point(0,0,0), normal(0,0,0), someface(-1) { }
	tTriangleMeshVertex(Point3 p) : point(p), normal(0,0,0), someface(-1) { }
	Point3 point;
	Vector3 normal;

	int someface;
};


class TriangleMeshFace {
public:

	TriangleMeshFace() {
		verts[0] = verts[1] = verts[2] = -1;
		nbrs[0]  = nbrs[1]  = nbrs[2]  = -1;
	}

	TriangleMeshFace(const int v[3]) {
		verts[0] = v[0];
		verts[1] = v[1];
		verts[2] = v[2];
		nbrs[0]  = nbrs[1] = nbrs[2] = -1;
	}

	int VertIndex(int v) const {
		if (verts[0] == v)	return 0;
		if (verts[1] == v)	return 1;
		if (verts[2] == v)	return 2;
		return -1;
	}

	int EdgeIndexCCW(int v0, int v1) const {

		int vi = VertIndex(v0);
		if (vi<0) return -1;

		if (verts[(vi+1)%3] == v1)
			return vi;

		return -1;
	}


	int verts[3];	// ccw order
	int nbrs[3];	// nbrs[0] shares verts 0,1, 1 shares 1,2, 2 shares 2,0
                    // -1 - boundary
};



template <class T>
class tTriangleMesh : public tModel<T>  {

public:
    typedef T value_type;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;
    typedef tBox3<T> Box3;
	typedef tTriangleMeshVertex<T> TriangleMeshVertex;

    typedef std::vector<TriangleMeshVertex> vertex_list;
    typedef std::vector<TriangleMeshFace> face_list;


	tTriangleMesh();
	~tTriangleMesh();

    bool Write(const char *fname) const;
    bool Read(const char *fname);

    bool ReadSMA(const char *fname);

    bool ReadOFF(const char *fname);
    bool WriteOFF(const char *fname) const; 

	// read a Hugues Hoppe .m format mesh file
	// we don't handle wedge/edge info
	bool ReadHHM(const char *fname, bool oside=false);
	bool WriteHHM(const char *fname) const;

    unsigned num_vertices() const;
    unsigned num_faces() const;
    TriangleMeshVertex& vertex(unsigned idx);
    const TriangleMeshVertex& vertex(unsigned idx) const;
    TriangleMeshFace& face(unsigned idx);
    const TriangleMeshFace& face(unsigned idx) const;

    bool VertOnBoundary(int v) const;
    void GetBoundaries(std::vector< std::vector<int> > &boundaries) const;

	// returns if the vert is on the boundary
	bool Vert1Ring(int v, std::vector<int> &ring) const;


	// do one iteration of loop subdivision
	void LoopSubdivide(bool keep_current_positions=false);

	void invalidate() {
	    tModel<T>::invalidate_all();
	}


#if 0
    //
    // Read/write a .off file
    //
    bool ReadOFF(const char* fname);
    void WriteOFF(const char* fname);
#endif

	// clear the whole mesh
	void Clear();

	void BuildStrips(int maxlen=0x7fffffff);


	// these are dangerous if the mesh has been changed!
	int FaceIndex(const TriangleMeshFace &f) const { return (&f - &faces[0]); }
	int VertexIndex(const TriangleMeshVertex &v) const { return (&v - &verts[0]); }

	Vector3 FaceNormal(const TriangleMeshFace &f) const {
		Vector3 e1 = verts[f.verts[1]].point - verts[f.verts[0]].point;
		Vector3 e2 = verts[f.verts[2]].point - verts[f.verts[0]].point;
		Vector3 norm = e1.cross(e2);
		norm.normalize();
		return norm;
	}

	void SetNormals();

	// Model stuff
	void compute_bounding_box() const;
	void compute_centroid() const;


	// iterator stuff
	class VertexVertexIterator;
	class VertexVertexIteratorI;

	class VertexFaceIterator {
	public:

		VertexFaceIterator(tTriangleMesh<T> &_mesh, int vi);
		bool done();
		VertexFaceIterator& operator++();
		TriangleMeshFace& operator*();

	private:

		bool first;
		int cur_face;
		int first_face;
		int vertex;
		tTriangleMesh<T> *mesh;

	    friend class tTriangleMesh<T>::VertexVertexIterator;
	};

	class VertexFaceIteratorI {
	public:

		VertexFaceIteratorI(const tTriangleMesh<T> &_mesh, int vi);
		bool done();
		VertexFaceIteratorI& operator++();
		int operator*();

	private:

		bool first;
		int cur_face;
		int first_face;
		int vertex;
		const tTriangleMesh<T> *mesh;

	    friend class tTriangleMesh<T>::VertexVertexIteratorI;
	};


	class VertexVertexIterator {
	public:

		VertexVertexIterator(tTriangleMesh<T> &_mesh, int vi);
		bool done();
		VertexVertexIterator& operator++();
		TriangleMeshVertex& operator*();
		bool isonedge() { return onedge; }

	private:

		bool first;
		bool onedge;
		VertexFaceIterator vfi;
	};


	class VertexVertexIteratorI {
	public:

		VertexVertexIteratorI(const tTriangleMesh<T> &_mesh, int vi);
		bool done();
		VertexVertexIteratorI& operator++();
		int operator*();
		bool isonedge() { return onedge; }

	private:

		bool first;
		bool onedge;
		VertexFaceIteratorI vfi;
	};


 //   void build_structures(std::map<int, int>& facemap,	std::map<int, int> vertmap);
	void build_structures(const std::vector<int> &facemap, const std::vector<int> &vertmap);
	void IdentityMap(std::vector<int> &map, int size);

//protected:
    vertex_list	verts;
    face_list	faces;

	std::vector<int> tstrips;	// triangle strip encoding for rendering


    bool read_a_line(FILE* f, char* line, int maxlinelen);

private:
	void BuildStrip(int f, std::vector<signed char> &face_avail, std::vector<int> &todo, int maxlen);

};

typedef tTriangleMesh<float> TriangleMeshf;
typedef tTriangleMesh<double> TriangleMeshd;

typedef tTriangleMeshVertex<float> TriangleMeshVertexf;
typedef tTriangleMeshVertex<double> TriangleMeshVertexd;

#if defined(REAL_IS_FLOAT)
typedef TriangleMeshf TriangleMesh;
#else
typedef TriangleMeshd TriangleMesh;
#endif

/*------------------------ Tools --------------------------*/

/*
 * For I/O, support conversion to and from IndexedTriangleSet 
 */
template <class T>
void TM2ITS(const tTriangleMesh<T>& tm, tIndexedTriangleSet<T>& its);
template <class T>
void ITS2TM(const tIndexedTriangleSet<T>& its, tTriangleMesh<T>& tm);

GTB_END_NAMESPACE

#include "triangle_mesh.inl"

#endif // __TRIANGLE_MESH_H
