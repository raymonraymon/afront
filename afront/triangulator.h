
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


#ifndef _TRIANGULATOR_H
#define _TRIANGULATOR_H

#include "front.h"
#include "uheap.h"
#include <gtb/gtb.hpp>
#include "guidance.h"


#define PROJECT_INCOMPLETE     -1       // projector thread hasn't gotten around to it yet
#define PROJECT_SUCCESS		0
#define PROJECT_FAILURE		1
#define PROJECT_BOUNDARY	2	// failed because we tried crossing a boundary


// this interface must be implemented and passed to the triangulator
// this stuff depends on the type of surface you are working on
class TriangulatorController {
    public:
    virtual ~TriangulatorController() {};

    // functions for storing the mesh as it's created
    // may want to save it for debugging drawing, or just write it out as it is created for more streamlike precessing
    virtual void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary) = 0;
    virtual void AddTriangle(int index, int v1, int v2, int v3) = 0;
    virtual void FinalizeVertex(int index) = 0;
    virtual void Finish() = 0;

    // how far are we allowed to step from this spot?
    virtual real_type MaxStepLength(const Point3 &p) = 0;
    virtual void InsertOverlayPoint(const Point3 &p, real_type curvature) = 0;
	virtual real_type LocalStepLength(real_type curvature) = 0;
    virtual void ResampleCurve(vector<Point3> &ip, vector< vector<Vector3> > &in,
							   vector<Point3> &op, vector< vector<Vector3> > &on,
							   bool is_not_loop) = 0;



    // project a point onto the surface
    virtual int ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const = 0;
    virtual int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const = 0;

    // get how big the fence should be to guarantee that there is no front crossing
    virtual real_type GetFenceScale() = 0;

    // get the next block to triangulate
    virtual void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const = 0;
    virtual bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const = 0;
	virtual int  GetBlockForPoint(const Point3 &p) const = 0;

};

// project a point onto the surface
class SurfaceProjector 
{
    public:
    virtual ~SurfaceProjector() {};
    virtual int ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const = 0;
    virtual int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const = 0;
    virtual void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const = 0;
    virtual bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const = 0;
	virtual int  GetBlockForPoint(const Point3 &p) const = 0;
};




////////////////////////////////////////////////////////////////////////////////

class OutputController
{
    public:
    OutputController() { child=NULL; }
    virtual ~OutputController() { if (child) delete child; }

    virtual void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary) = 0;
    virtual void AddTriangle(int index, int v1, int v2, int v3) = 0;
    virtual void FinalizeVertex(int index) = 0;
    virtual void Finish() = 0;


    static void AddControllerToBack(OutputController *&h, OutputController *c) {
		if (!h) {
			h = c;
		} else {
			AddControllerToBack(h->child, c);
		}
    }
    OutputController *child;
};



class OutputControllerITS : public OutputController
{
    public:
    OutputControllerITS() { }
    virtual ~OutputControllerITS() { };

    void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary) {
		triangulation.add_vertex(p);
		triangulation.insert_vertex_normal(n);

		if (child)
			child->AddVertex(index, p, n, boundary);
    }
  
    void AddTriangle(int index, int v1, int v2, int v3) {
		triangulation.insert_triangle((unsigned)v1, (unsigned)v2, (unsigned)v3, false, 0.0);

		if (child)
			child->AddTriangle(index, v1, v2, v3);
    }

    void FinalizeVertex(int index) {
		if (child)
			child->FinalizeVertex(index);
    }


    void Finish() {
		if (child)
			child->Finish();
    }

    IndexedTriangleSet triangulation;
};


////////////////////////////////////////////////////////////////////////////////
// wrapper to allow different implementations of output, guidance field, and surface

class ControllerWrapper : public TriangulatorController 
{
    public:
    ControllerWrapper(GuidanceField *g, SurfaceProjector *p, OutputController *oc)
		: TriangulatorController(),
		  guidance(g), projector(p), output_c(oc) { }

    virtual ~ControllerWrapper() {
    }

    void SetProjector(SurfaceProjector *p) {
		projector = p;
    }


    void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary) {
		output_c->AddVertex(index, p, n, boundary);
    }

    void AddTriangle(int index, int v1, int v2, int v3) {
		output_c->AddTriangle(index, v1, v2, v3);
    }

    void FinalizeVertex(int index) {
		output_c->FinalizeVertex(index);
    }

    void Finish() {
		output_c->Finish();
    }


    real_type MaxStepLength(const Point3 &p) {
		return guidance->MaxStepLength(p);
    }

    void InsertOverlayPoint(const Point3 &p, real_type curvature) {
		guidance->InsertOverlayPoint(p, curvature);
	}


	real_type LocalStepLength(real_type curvature) {
		return guidance->LocalStepLength(curvature);
	}


    void ResampleCurve(vector<Point3> &ip, vector< vector<Vector3> > &in,
					   vector<Point3> &op, vector< vector<Vector3> > &on,
					   bool is_not_loop) {
		guidance->ResampleCurve(ip, in, op, on, is_not_loop);
    }



    int ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const {
		return projector->ProjectPoint(base1, base2, fp, fn, tp, tn, curvature);
    }

    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const {
		return projector->ProjectPoint(fp, tp, tn);
    }

    real_type GetFenceScale() { return guidance->GetFenceScale(); }


    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
		projector->GetFirstBlock(box, pts, norms, loops, seeds);
    }

    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
		return projector->GetNextBlock(box, pts, norms, loops, seeds);
    }

	int GetBlockForPoint(const Point3 &p) const {
		return projector->GetBlockForPoint(p);
	}

    private:

    GuidanceField *guidance;
    SurfaceProjector *projector;
    OutputController* output_c;
};



class Triangulator {

    public:
    // stuff needed by the kd tree
    typedef ::Box3 Box3;
    typedef ::Point3 Point3;


    typedef Front::feli feli;
    typedef gtb::BoxKDTree<feli, Triangulator> kdtree_type;
    typedef UHeap<feli> heap_type;


    Triangulator(TriangulatorController &c);
    ~Triangulator();


    void InsertSubMesh(const TriangleMesh &mesh, const vector<int> &sides, int sidetoadd);
    /*!
     * Triangulate a region
     * Input:
     *    ipts - set of points that form the initial front?
     *    inorms - their normals.
     */
    void Go(const vector< vector<Point3> > &ipts, const vector< vector<Vector3> > &inorms, 
			bool failsafe,
			const vector<Point3> &crease_points,
			const vector<Vector3> &crease_onormals,
			const vector< vector<int> > &crease_indices,
			const vector< vector<Vector3> > &crease_normals,
			const vector< int > &corner_triangles);
    void Go(const vector< vector<Point3> > &ipts, const vector< vector<Vector3> > &inorms, bool dofailsafe);
    void ExtractFronts(vector< vector<Point3> > &fpts,
					   vector< vector<Vector3> > &fnorms,
					   vector< vector<int> > &states,
					   vector< vector<real_type> > &fenceheight);


    Box3 bounding_box(const feli &i) const;
    real_type distance(const feli &i, const Point3 &from) const;

    void SetFlipOutput(bool f) { flipOutput=f; }

	void SetOutputIndices(int v, int f) { numVertsAdded=v; numFacesAdded=f; }
	void GetOutputIndices(int &v, int &f) const { v=numVertsAdded; f=numFacesAdded; }



    private:

    void CreateVertex(const Point3 &p, const Vector3 &n, FrontElement &fe, bool boundary, real_type step);
    void CreateTriangle(const feli v1, const feli v2, const feli v3, bool forreal=true);
    void CreateTriangle(int v1, int v2, int v3);

    void PrioritizeEdgeGrow(feli e1);
    void PrioritizeEdgeConnect(feli e1);
    void PrioritizeEdgeFailsafe(feli e);
    bool EdgeExistsAlready(feli e1, feli e2) const;

    bool GetConnection(feli e1, feli &across);
    real_type CosSmallestRemainingAngle(feli e1, feli e2, feli across);


    real_type PrioritizeConnection(const FrontElement &p1, const FrontElement &p2, const FrontElement &p3);
    real_type PrioritizeEdgeGrow(const FrontElement &p1, const FrontElement &p2);


    bool IntersectsFence(const feli &e, const Point3 &p1, const Point3 &p2) const;
    bool InLegalFenceCorner(const feli &e, const Point3 &p) const;
    real_type DistanceToFence(const feli &e1, const Point3 &p) const;
    real_type DistanceToFence(const Point3 &ep1, const Point3 &ep2, const Vector3 &en1, const Vector3 &en2, const Point3 &p) const;
    bool GrowthFromEdgeLegal(const feli e1, const Point3 &p) const;
    bool TriangleLegal(const feli e1, const feli e2, const feli *across, const Point3 &across_p, const Vector3 &across_n, vector<feli> *possible_intersects, bool *isclose) const;

    void ConnectTriangle(feli e1, feli across, bool dofailsafe);
    void GrowEdge(feli e1, const Point3 &p, const Vector3 &v, real_type step);


    class ProjectionWork {
		public:
		ProjectionWork() { };
		ProjectionWork(real_type _priority, ProjectionResult *_pr, FrontElement *_base1, const FrontElement *_base2) :
			priority(_priority), pr(_pr), base1(_base1), base2(_base2) { }
		bool operator<(const ProjectionWork &rhs) const {
			return (priority < rhs.priority);
		}
		real_type priority;
		ProjectionResult *pr;
		const FrontElement *base1;
		const FrontElement *base2;
    };
    void RequestProjection(feli e);
    void WaitForProjection(feli e);
    void CancelProjection(feli e);
    static void* ProjectorThreadMain(void *arg);
    void PauseProjections();
    void ContinueProjections();
	void DoProjection(ProjectionWork &pw);
    static bool GetTentativePoint(const Point3 &p1, const Vector3 &n1, real_type s1,
								  const Point3 &p2, const Vector3 &n2, real_type s2,
								  Point3 &tentative_p, Vector3 &tentative_n);


    gtb::fast_pq<ProjectionWork> work_queue;
    volatile bool work_pause;
    volatile int work_num_paused;
    volatile bool work_quit;
    CSObject work_cs;
    vector<thlib::Thread*> work_threads;
	

    void StartWorkerThreads();
    void StopWorkerThreads();


    void SetupInitialCreaseFronts
		(const vector<Point3> &points, // The points on the creases
		 const vector<Vector3> &onormals, // The normals for the output
		 const vector< vector<int> > &indices, // The front element point indices
		 const vector< vector<Vector3> > &normals, // The normals for the front elements
		 const vector< int > &corner_triangles); 

    void SetupInitialFront(const vector<Point3> &ipts, const vector<Vector3> &inorms, bool loop);
    void SetupInitialFronts(const vector< vector<Point3> > &ipts, const vector< vector<Vector3> > &inorms);
    void RemoveFrontElement(feli e);

    bool InWorkingArea(feli e1, feli e2);
    void UpdateWorkingArea(feli e);
    real_type workingAreaRadius;
    Point3 workingAreaCenter;

    bool InWorkingBlock(feli e1, feli e2);
    void UpdateWorkingBlock(bool first=false);
    bool haveWorkingBlock;
    Box3 workingBlock;

    void InsertBoundary(vector<Point3> &pts, vector<Vector3> &norms, bool isloop);
    vector<feli> boundary_connectors;

    void VerifyFronts();

    void VertexRefIncrement(FrontElement &v);
    void VertexRefDecrement(FrontElement &v);

    void AddNeighborsToRings(FrontElement &v1, FrontElement &v2);

    class VertInfo {
		public:
		int front_refcount;
		int vertindex;
		vector<int> ring;
    };

    vector<VertInfo> vertinfo;
    vector<int> available_vertinfo;
    int AllocVertInfo();
    void FreeVertInfo(int i);
    int max_active_verts;
    

	void KDInsert(feli e);
	void KDRemove(feli e);
	void KDIntersectedBoxes(const Box3 &box, vector<feli> &intersects) const;
	void KDStartOT(const Point3 &p);
	real_type KDNextOT(feli &next);
	void KDEndOT();
	void KDUpdate(feli e);
	
	feli next_l, next_g;
	real_type dnext_l, dnext_g;
	kdtree_type::OrderedTraverse *ot_l, *ot_g;


    kdtree_type kdtree_l;
    kdtree_type kdtree_g;

    heap_type	heap;
    vector<Front*> fronts;

    TriangulatorController &controller;

    int numVertsAdded;
    int numFacesAdded;
    bool flipOutput;
	
};





#endif

