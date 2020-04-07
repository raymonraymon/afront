
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


// why aren't these dependent on the reduction factor?
static const float cos_ear_cuttable_angle_max = cosf(M_PI * 100.0f / 180.0f);
static const float cos_ear_cuttable_angle_min = cosf(M_PI * 0.0f / 180.0f);





///////////////////////////////////////////////////////////////////////////////
// GuidanceField


real_type GuidanceField::StepRequired(real_type dist, real_type ideal) const
{
    return (1-reduction)*dist + reduction*ideal;
}


real_type GuidanceField::MaxCurvatureToIdeal(real_type curvature) const {
    return 2*sin(rho*0.5) / fabs(curvature);
}


real_type GuidanceField::LocalStepLength(real_type curvature) {
	return StepRequired(0, MaxCurvatureToIdeal(curvature));
}


real_type GuidanceField::MaxStepLength(const Point3 &p)
{
    real_type len = max_step;
    real_type checked_rad=0;

    void *ctx = OrderedPointTraverseStart(p);
    do {
		
		Point3 point;
		real_type ideal;
		if (!OrderedPointTraverseNext(ctx, checked_rad, point, ideal)) break;


		// apparently the distance returned by the ordered traverse isn't the actual distance!
		checked_rad = sqrt(checked_rad);
		len = std::min(len, StepRequired(checked_rad, ideal));
		if (len<min_step) {
			len = min_step;
			break;
		}

    } while (checked_rad < (len/(1-reduction)));
    OrderedPointTraverseEnd(ctx);

	// check the overlay
	// only if we really need to - it has a synchronization lock in it!
	extern bool retro;
	if (retro)
		len = o_MaxStepLength(p, len);

    return len;
}


real_type GuidanceField::o_MaxStepLength(const Point3 &p, real_type startlen)
{
    real_type len = startlen;
    real_type checked_rad=0;

    void *ctx = o_OrderedPointTraverseStart(p);
    do {
		
		Point3 point;
		real_type ideal;
		if (!o_OrderedPointTraverseNext(ctx, checked_rad, point, ideal)) break;

		// apparently the distance returned by the ordered traverse isn't the actual distance!
		checked_rad = sqrt(checked_rad);
		len = std::min(len, StepRequired(checked_rad, ideal));
		if (len<min_step) {
			len = min_step;
			break;
		}

    } while (checked_rad < (len/(1-reduction)));


	o_OrderedPointTraverseEnd(ctx);

    return len;
}


void GuidanceField::ResampleCurve(const vector<Point3> &ip, const vector< vector<Vector3> > &in,
								  vector<Point3> &op, vector< vector<Vector3> > &on,
								  bool is_not_loop) {

    real_type totallength=0;
    real_type steplengths=0;
    for (unsigned i=0; i<ip.size() - (int)is_not_loop; i++) {
		totallength += Point3::distance(ip[i], ip[(i+1)%ip.size()]);
    }

    vector<real_type> steps;

    for (int pass=0; pass<2; pass++) {

		op.clear();
		op.push_back(ip[0]);

		on.resize(in.size());
		for (unsigned i=0; i<in.size(); i++) {
			on[i].clear();
			on[i].push_back(in[i][0]);
		}


		int lastpassed = 0;
		bool done = false;
		while (!done) {

			if (pass==1 && op.size()==steps.size())
				break;

			Point3 lp = op.back();

			real_type dist_to_go;
			if (pass==0) {
				steps.push_back(MaxStepLength(lp));
				steplengths += steps[op.size()-1];
				dist_to_go = steps[op.size()-1];
			} else {
				dist_to_go = steps[op.size()-1] * totallength / steplengths;
				//				dist_to_go = steps[op.size()-1] * steplengths / (steplengths-steps.back());
			}


			while (dist_to_go > 0) {

				int nextpoint = (lastpassed+1) % ip.size();

				real_type pdist = Point3::distance(lp, ip[nextpoint]);

				if (pdist >= dist_to_go) {

					// don't go all the way to the endpoint

					Vector3 edir = ip[nextpoint] - ip[(nextpoint+ip.size()-1)%ip.size()];
					Vector3 edirn = edir.normalized();
					op.push_back(lp + dist_to_go * edirn);

					real_type alpha = Point3::distance(op.back(), ip[nextpoint]) / edir.length();

					for (unsigned i=0; i<in.size(); i++) {
						on[i].push_back(alpha*in[i][(nextpoint+ip.size()-1)%ip.size()] + (1-alpha)*in[i][nextpoint]);
						on[i].back().normalize();
					}

					dist_to_go = 0;

				} else {

					dist_to_go -= pdist;
					lp = ip[nextpoint];
					lastpassed = nextpoint;
					if ((is_not_loop && lastpassed==ip.size()-1) ||
						(!is_not_loop && lastpassed == 0)) {
						dist_to_go = 0;
						done = true;
					}
				}
			}
		}
		if (is_not_loop) {
			op.push_back(ip.back());
			for (unsigned i=0; i<in.size(); i++) {
				on[i].push_back(in[i].back());
			}
		}
    }
}



void* GuidanceField::o_OrderedPointTraverseStart(const Point3 &p) {

	while (1) {
		o_sl.enter();

		if (!o_pause) {
			o_nopt++;
			o_sl.leave();
			break;
		}
		
		o_sl.leave();
		//		thlib::Thread::yield();
	}


    if (o_kdtree) {
		return new o_kdtree_type::OrderedIncrementalTraverse(*o_kdtree, p);
	} else {
		return NULL;
	}
}


bool GuidanceField::o_OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
    if (!ctx) {
		return false;
	}

	o_kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (o_kdtree_type::OrderedIncrementalTraverse*)ctx;
	if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = o_gfpoints[n];
    ideal = o_ideal_lengths[n];
    return true;
}


void GuidanceField::o_OrderedPointTraverseEnd(void *ctx) {
	if (ctx) {
		o_kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (o_kdtree_type::OrderedIncrementalTraverse*)ctx;
		delete kdOrderedTraverse;
	}

	// decrement the nopt count
	o_sl.enter();
	o_nopt--;
	o_sl.leave();
}


void GuidanceField::InsertOverlayPoint(const Point3 &p, real_type curvature) {

	extern bool retro;
	if (!retro) {
		cerr<<"inserting overlay point without retro!!"<<endl;
		exit(1);
	}

	o_sl.enter();

	// check if someone else is already trying to add a point 
	while (o_pause) {
		o_sl.leave();
		//		thlib::Thread::yield();
		o_sl.enter();
	}

	o_pause = true;

	// wait for any queries to finish
	while (o_nopt != 0) {
		o_sl.leave();
		//		thlib::Thread::yield();
		o_sl.enter();
	}


	o_gfpoints.push_back(p);
	o_ideal_lengths.push_back(MaxCurvatureToIdeal(curvature));

	// rebuild the tree if we need to
	if (o_gfpoints.size() > o_lastkdsize*2) {
		cerr<<"rebuilding overlay tree with "<<o_gfpoints.size()<<" points"<<endl;

		if (o_kdtree)
			delete o_kdtree;

		o_kdtree = new o_kdtree_type(10, Box3::bounding_box(o_gfpoints), o_kdGetPoint);
		for (unsigned i=0; i<o_gfpoints.size(); i++)
			o_kdtree->Insert(i);
		o_kdtree->MakeTree();

		o_lastkdsize = o_gfpoints.size();

		cerr<<"ok"<<endl;
	} else {
		// just insert
		o_kdtree->Insert(o_gfpoints.size()-1);
	}


	// let things continue
	o_pause = false;
	o_sl.leave();
}


