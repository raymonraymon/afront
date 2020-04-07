
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



#ifndef _TRIANGULATE_CSG
#define	_TRIANGULATE_CSG


// for the pointsides arrays:
// -2: point is close to the intersection on the side to keep
// -1: point is far from the intersection on the side to keep
//  1: point is far from the intersection on the side to get rid of
//  2: point is close to the intersection on the side to get rid of




class TriangleMeshFaceTree {
    public:

    typedef ::Box3 Box3;
    typedef ::Point3 Point3;

    TriangleMeshFaceTree(const TriangleMesh &m) : mesh(m) { }
    Box3 bounding_box(int i) const {
	return Box3::bounding_box(mesh.verts[mesh.faces[i].verts[0]].point,
				  mesh.verts[mesh.faces[i].verts[1]].point,
				  mesh.verts[mesh.faces[i].verts[2]].point);
    }

    real_type distance(int i, const Point3 &p) const {
	Point3 closest = gtb::tTriangle3<real_type>(mesh.verts[mesh.faces[i].verts[0]].point,
						    mesh.verts[mesh.faces[i].verts[1]].point,
						    mesh.verts[mesh.faces[i].verts[2]].point).closest_point(p);
	return Point3::distance(p, closest);
    }

    const TriangleMesh &mesh;
};



class MeshCSGGuidanceField : public GuidanceField {
    public:

    MeshCSGGuidanceField(int curv_sub, const TriangleMesh &mesh1, const TriangleMesh &mesh2, const vector<int> pointsides[2], vector< vector<Point3> > &curves, real_type rho, real_type min_step, real_type max_step, real_type reduction);
    ~MeshCSGGuidanceField();

    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);

    void Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad);

    static void GetIntersectionLoops(const TriangleMesh &m1, const TriangleMesh &m2, vector< vector<Point3> > &loops, vector< vector<Vector3> > &n1, vector< vector<Vector3> > &n2, vector<int> pointsides[2]);
    static void FloodPointSides(const TriangleMesh &m, vector<int> &sides);

    void TrimPointSides(const TriangleMesh &mesh, vector<int> &sides);
    static void FixPointSides(const TriangleMesh &mesh, vector<int> &sides);
    static void FrontsFromSides(const TriangleMesh &mesh, vector<int> &sides, vector< vector<Point3> > &fpoints, vector< vector<Vector3> > &fnormals);
    void BlendGuidance(bool blend, const TriangleMesh* meshes[2], vector<int> sides[2]);



    private:
    class GetPoint {
	public:
	vector<Point3> allpoints;
	GetPoint() {}
	const Point3& operator()(unsigned idx) const {
	    return allpoints[idx];
	}
    };

    typedef gtb::KDTree<int, real_type, GetPoint> kdtree_type;
    kdtree_type *kdtree;
    GetPoint kdGetPoint;

    vector<real_type> ideal_length;

};



// for mesh/pointset mixed mode csg
class MeshPSCSGGuidanceField : public GuidanceField {
    public:

    MeshPSCSGGuidanceField(int curv_sub, const TriangleMesh &mesh, const surfel_set &points, CProjection &psprojector, vector< vector<Point3> > &curves, real_type rho, real_type min_step, real_type max_step, real_type reduction, int adamson);
    ~MeshPSCSGGuidanceField();

    void* OrderedPointTraverseStart(const Point3 &p);
    int OrderedPointTraverseNext(void *ctx, real_type &squared_dist);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);

    void Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad);


    static void GetIntersectionLoops(MeshProjector &mp, SmoothMLSProjector &psp,
				     vector< vector<Point3> > &loops, vector< vector<Vector3> > &mnorms, vector< vector<Vector3> > &pnorms);

    void InitPointSides(const TriangleMesh &m, const surfel_set &points, vector<int> &sides, const vector< vector<Point3> > &loops);

    static void HackPointSides(const TriangleMesh &m, vector<int> &sides);


    private:


    static void ProjectPointOntoIntersection(MeshProjector &mp, SmoothMLSProjector &psp,
					     const Point3 &from, const Point3 *last, real_type dist, Point3 &to, Vector3 &mnorm, Vector3 &pnorm);


    class GetPoint {
	public:
	vector<Point3> allpoints;
	GetPoint() {}
	const Point3& operator()(unsigned idx) const {
	    return allpoints[idx];
	}
    };

    typedef gtb::KDTree<int, real_type, GetPoint> kdtree_type;
    kdtree_type *kdtree;
    GetPoint kdGetPoint;

    vector<real_type> ideal_length;
};





#endif

