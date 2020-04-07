
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


#include "common.h"
#include "parallel.h"
#include "guidance.h"
#include "triangulator.h"
#include "triangulate_mls.h"

#include <cstdio>
#include <algorithm>


#define NO_RMLS
#ifndef NO_RMLS
#include <rmlslib/Primitive_Functions.h>
#endif

using namespace std;

SmoothMLSProjector::SmoothMLSProjector(surfel_set &surfels, int adamson):
    _wf(1.0f),
    _radius_wf(1.0f),
    _projector(surfels, 8, &_wf, &_radius_wf, 2, NULL),
    _adamson(adamson)
{
    _projector.compute_points_radius();
    _projector.set_radius_factor(2.0f);
}

SmoothMLSProjector::~SmoothMLSProjector()
{
}

int SmoothMLSProjector::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const
{
	return (_projector.PowellProject(fp, tp, tn) ? PROJECT_SUCCESS : PROJECT_FAILURE);
}

int SmoothMLSProjector::ProjectPoint
    (const FrontElement &base1, const FrontElement &base2, 
     const Point3 &fp, const Vector3 &fn,
     Point3 &tp, Vector3 &tn, real_type *curvature) const
{
    extern bool rmls;
    if (rmls) {
#ifndef NO_RMLS

	real_type radius = _projector.point_radius(fp);
	gtb::surfelset_view nbhd(_projector.get_points());

	extern int rmls_knn;
	_projector.extract2(fp, rmls_knn, nbhd);

	// copy into a vector for rmls
	vector<Point> neighbors;
	for (unsigned i=0; i<nbhd.size(); i++) {
	    const Point3 &p = nbhd.vertex(i);
	    neighbors.push_back(Point(p[0], p[1], p[2]));
	}

	Point rmls_fp(fp[0], fp[1], fp[2]);
	extern real_type noise_threshold;
	Vertex res = rmls_projection(rmls_fp, neighbors, noise_threshold);
	tp = Point3(res.point_[0], res.point_[1], res.point_[2]);
	tn = Vector3(res.normal_.x(), res.normal_.y(), res.normal_.z());

	// FIXME - is it possible to fail??
#else
		cerr<<"robust mls code not released yet!"<<endl;
		exit(0);
#endif
    } else {
		if (_adamson != 0) {
			if (!_projector.adamson_projection(fp, tp, tn, _adamson))
				return PROJECT_FAILURE;
		} else {
			if (!_projector.PowellProject(fp, tp, tn))
				return PROJECT_FAILURE;
		}
	}

    if (fn.dot(tn) < 0)
		tn.flip();
    if (fn.dot(tn) < 0.7) {
		//		cerr << "dUuuuuh" << endl;
		tn = fn;
    }
    return PROJECT_SUCCESS;
}


void SmoothMLSGuidanceField::curvatures(const Point3 &ref, Vector3 &normal, real_type &k1, real_type &k2) const
{
    Point3 result;
    real_type radius = _projector.point_radius(ref);
    surfelset_view nbhd(_projector.get_points());
    _projector.extract2(ref, radius, nbhd);

    if (_adamson != 0) {
	_projector.adamson_projection(nbhd, ref, result, normal, _adamson, &k1, &k2);
    } else {
	mls::Poly2<real_type> poly;
	plane_transformation T;
	surfel_set std_points;
	if (_projector.PowellProject(nbhd, ref, result, normal, T, &poly, std_points))
	    poly.curvatures(0, 0, k1, k2);
	else
	    k1=k2=1e-5;
    }

    //	cerr<<k1<<" "<<k2<<endl;
}

void SmoothMLSGuidanceField::MakeGuidanceParallel(int nt, int id, CSObject &cs, vector<real_type> &ideals) const
{
	for (unsigned i=id; i<_projector.get_points().size(); i+=nt) {
		real_type k1, k2;
		Vector3 normal;
		Point3 p = _projector.get_points().vertex(i);
		curvatures(p, normal, k1, k2);

		ideals[i] = MaxCurvatureToIdeal(k2);
			
		if ((int)(i * 100.0 / _projector.get_points().size()) != (int)((i+1) * 100.0 / _projector.get_points().size())) {
			//			cerr<<"\r          \r"<<((int)((i+1) * 100.0 / _projector.get_points().size()))<<"%";
			cerr<<".";
			cerr.flush();
		}
	}
}

SmoothMLSGuidanceField::SmoothMLSGuidanceField
    (CProjection &projector,
     real_type rho, real_type min_step, real_type max_step, real_type reduction, 
     int adamson, const char *gffile) :
	GuidanceField(rho, min_step, max_step, reduction),
	_projector(projector),
	_adamson(adamson)
{
    ideal_length.resize(projector.get_points().size());
    cerr << "Point count: ";
    cerr << projector.get_points().size() << endl;


	if (gffile != NULL) {
		cerr << "reading guidance field from file " << gffile;
		
		FILE *gff = fopen(gffile, "r");
		if (gff == NULL) {
			cerr << "couldn't open guidance field file " << gffile;
			exit(0);
		}

		for (unsigned i=0; i<projector.get_points().size(); i++) {
			char line[1000];
			fgets(line, 1000, gff);

			float ideal = 0;
			sscanf(line, "%g", &ideal);
			ideal_length[i] = atof(line);
		}

		fclose(gff);
	}

	else {
		ParallelExecutor(idealNumThreads, makeClassFunctor(this,&SmoothMLSGuidanceField::MakeGuidanceParallel), ideal_length);


		extern int mls_gf_fix;
		for (int it=0; it<mls_gf_fix; it++) {

			for (int i=0; i<_projector.get_points().size(); i++) {
				surfelset_view nbhd(_projector.get_points());
				_projector.extract(_projector.get_points().vertex(i),
								   ((real_type)1)*_projector.point_radius(_projector.get_points().vertex(i)),
								   nbhd);

				real_type smallest = INFINITY;
				for (int j=0; j<nbhd.size(); j++) {
					if (fabs(ideal_length[nbhd.get_index(j)]) < fabs(smallest))
						smallest = ideal_length[nbhd.get_index(j)];
				}

				if (smallest == ideal_length[i])
					ideal_length[i] = INFINITY;
			}			
		}
	}

	cerr << endl;
}

SmoothMLSGuidanceField::~SmoothMLSGuidanceField()
{
}


void* SmoothMLSGuidanceField::OrderedPointTraverseStart(const Point3 &p) {
    return new gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse(*_projector.get_kdtree().tree, p);
}


bool SmoothMLSGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
	gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse *kdOrderedTraverse = (gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse*)ctx;

    if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = _projector.get_points().vertex(n);
    ideal = ideal_length[n];
    return true;
}


void SmoothMLSGuidanceField::OrderedPointTraverseEnd(void *ctx) {
	gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse *kdOrderedTraverse = (gtb::ss_kdtree<surfel_set>::t_surfel_tree::OrderedIncrementalTraverse*)ctx;
    delete kdOrderedTraverse;
}


void SmoothMLSGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad)
{
    pts = _projector.get_points().vertices();
    norms = _projector.get_points().normals();
    rad = ideal_length;
}
