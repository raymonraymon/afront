
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



#ifndef _TRIANGULATE_MESH
#define	_TRIANGULATE_MESH

#include <viewer/debug_draw.h>
#include "common.h"
#include "guidance.h"
#include "triangulator.h"

class MeshProjector : public SurfaceProjector {
    public:
    // stuff for the kdtree
    typedef ::Box3 Box3;
    typedef ::Point3 Point3;

    MeshProjector(const TriangleMesh &mesh, const vector<int> *pointsides=NULL);
    virtual ~MeshProjector();
    int ProjectPoint(const FrontElement &base1, const FrontElement &base2, 
		     const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const;

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
	box = mesh.bounding_box();
	box.enlarge(0.1);
    }
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const { return false; }
	int GetBlockForPoint(const Point3 &p) const { return 0; }

    Box3 bounding_box(int face) const;
    real_type distance(int face, const Point3 &from) const;


    // this is a simple closest point projection
    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const;
    typedef gtb::BoxKDTree<int, MeshProjector> kdtree_type;
    kdtree_type *GetKdTree();

    
    void FindSoupBoundaries(vector< vector<int> > &boundaries, real_type d) const;

    private:

    const TriangleMesh &mesh;
    kdtree_type *kdtree;
};





class MeshGuidanceField : public GuidanceField {
    public:

    MeshGuidanceField(int curv_sub, const TriangleMesh &mesh, real_type rho, real_type min_step, real_type max_step, real_type reduction);
    virtual ~MeshGuidanceField();

    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);

    real_type IdealLength(int i) const { return ideal_length[i]; }

    void Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad);

    class GetPointMesh {
	public:
	const TriangleMesh& mesh;
	GetPointMesh(const TriangleMesh& m) : mesh(m) {}
	const Point3& operator()(unsigned idx) const {
	    return mesh.verts[idx].point;
	}
    };
    typedef gtb::KDTree<int, real_type, GetPointMesh> kdtree_type;
    kdtree_type *GetKdTree();


    private:

    vector<real_type> ideal_length;

    kdtree_type *kdtree;
    GetPointMesh kdGetPoint;
};



// compute curvature of a space curve by fitting a quadratic polynomial
real_type edge_curvature(const vector<Point3> &points, const vector<real_type> &weights, int si);



#endif

