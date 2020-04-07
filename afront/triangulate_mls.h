
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


#ifndef _TRIANGULATE_MLS_SMOOTH
#define _TRIANGULATE_MLS_SMOOTH

#include <gtb/gtb.hpp>
#include <mlslib/CProjection.h>
#include <vector>
#include "common.h"

class SmoothMLSProjector : public SurfaceProjector 
{
    public:
    SmoothMLSProjector(surfel_set &surfels, int adamson);
    virtual ~SmoothMLSProjector();
    
    mls::GaussianWeightFunction<real_type> _wf, _radius_wf;
    mls::CProjection<real_type> _projector;
    int _adamson;

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
	box = _projector.get_points().bounding_box();
	box.enlarge(0.1);
    }
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const { return false; }
	int GetBlockForPoint(const Point3 &p) const { return 0; }


    int ProjectPoint(const FrontElement &base1, const FrontElement &base2,
		     const Point3 &fp, const Vector3 &fn,
		     Point3 &tp, Vector3 &tn, real_type *curvature) const;

    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const;
    CProjection &get_projector();

    void SaveField(const char *filename);
    void SetRadiusFactor(real_type t) {
	_projector.set_radius_factor(t);
    }
};

class SmoothMLSGuidanceField : public GuidanceField 
{
    public:
    typedef Vector3 (*VectorField)(const Point3 &);

    // Compute from guidance field
    SmoothMLSGuidanceField(CProjection &projector, 
						   real_type rho,
						   real_type min_step, real_type max_step, real_type reduction, int adamson,
						   const char *gffile=NULL);
    virtual ~SmoothMLSGuidanceField();

	void MakeGuidanceParallel(int nt, int id, CSObject &cs, vector<real_type> &ideals) const;


    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);

    real_type IdealLength(int i) const { return ideal_length[i]; }

    void Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad);

    void curvatures(const Point3 &ref, Vector3 &normal, real_type &k1, real_type &k2) const;

    CProjection &_projector;

    gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse *_kdOrderedTraverse;

    int _adamson;

    private:

    vector<real_type> ideal_length;

};

#endif
