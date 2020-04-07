
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



#ifndef _GUIDANCE_H
#define _GUIDANCE_H

#include "front.h"

#include <gtb/graphics/point3.hpp>

#include "conekdtree.h"


template <typename T>
	class GetPointVector {
    public:
    GetPointVector(const vector<T> &pts) : _pts(pts) { }
    const T& operator()(unsigned idx) const { return _pts[idx]; }
    const vector<T> &_pts;
};



class GuidanceField {
    public:

    GuidanceField(real_type _rho, real_type _min, real_type _max, real_type _reduction)
		: rho(_rho), min_step(_min), max_step(_max), reduction(_reduction),
		  o_kdGetPoint(o_gfpoints), o_kdtree(NULL), o_nopt(0), o_pause(false), o_lastkdsize(0) { };
    virtual ~GuidanceField() {};


    virtual void* OrderedPointTraverseStart(const Point3 &p) = 0;
    virtual bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) = 0;
    virtual void OrderedPointTraverseEnd(void *ctx) = 0;

    real_type MaxStepLength(const Point3 &p);
	real_type LocalStepLength(real_type curvature);
	void InsertOverlayPoint(const Point3 &p, real_type curvature);

    void ResampleCurve(const vector<Point3> &ip, const vector< vector<Vector3> > &in,
					   vector<Point3> &op, vector< vector<Vector3> > &on,
					   bool is_not_loop=false);

    real_type GetFenceScale() {
		extern real_type fence_scale;
		return (real_type)(fence_scale * 4.0*((1-sqrt((1+2*cos(rho))/3)) / (2*sin(rho*0.5))));
    }


    real_type MaxCurvatureToIdeal(real_type curvature) const;


    template <typename GETPOINT>
		class Trimmer {
		public:
		Trimmer(const GuidanceField &g, const GETPOINT &gp, const vector<real_type> &id) :
			guidance(g), getpoint(gp), ideal(id) { }


		void Trim(vector<int> &marked) const {

			if (marked.size() == 0)
				return;

			for (unsigned i=0; i<marked.size(); i++) {
				marked[i] = 0;
			}

			extern bool trim_guidance;
			if (!trim_guidance) {
				return; 
			}

			vector<int> ipts(marked.size());
			for (unsigned i=0; i<ipts.size(); i++) {
				ipts[i] = i;
			}


			RecursiveTrim(ipts, marked);
			int nmarked = marked.size() - std::accumulate(marked.begin(), marked.end(), 0);

			ipts.resize(0);
			for (unsigned i=0; i<marked.size(); i++) {
				if (!marked[i])
					ipts.push_back(i);
			}

			RecursiveTrim(ipts, marked);
			nmarked = marked.size() - std::accumulate(marked.begin(), marked.end(), 0);
		}


		// stuff for the cone kd-tree
		Box4 bounding_box(int v) const {
			real_type r = guidance.StepRequired(0, ideal[v]);
			Point3 p = getpoint(v);
			return Box4(Box3(p,p), r, r);
		}



		private:

		void ParallelTrim(int nt, int id, CSObject &cs,
						  vector< std::pair<real_type,int> > &points,
						  ConeBoxKDTree<int, Trimmer> &kd,
						  vector<int> &marked) const {
			real_type t = 1-guidance.reduction;
			t *= t;
			real_type angle_dot = sqrt(t / (t+1));

			unsigned count = points.size();
			for (unsigned i=id; i<count; i+=nt) {
				if (marked[points[i].second])
					continue;
				Cone c;
				c.p = getpoint(points[i].second);
				c.r = guidance.StepRequired(0,ideal[points[i].second]);
				c.angle_dot = angle_dot;
				kd.MarkIntersected(*this, c, marked);
			}
		}


		void RecursiveTrim(const vector<int> &ipts, vector<int> &marked) const {

			extern int trim_bin_size;
			if ((int)ipts.size() < trim_bin_size) {

				vector< std::pair<real_type,int> >  sorted(ipts.size());
				for (unsigned i=0; i<ipts.size(); ++i) {
					sorted[i].first = guidance.StepRequired(0,ideal[ipts[i]]);
					sorted[i].second = ipts[i];
				}
				sort(sorted.begin(), sorted.end());

				ConeBoxKDTree<int, Trimmer> kd(ipts, *this);
				ParallelExecutor(idealNumThreads, makeClassFunctor(this, &Trimmer::ParallelTrim), sorted, kd, marked);

			} else {
				// compute the bounding box
				Box3 bbox(getpoint(ipts[0]), getpoint(ipts[0]));
				for (unsigned i=1; i<ipts.size(); i++) {
					bbox.update(getpoint(ipts[i]));
				}

				vector<int> sides[2];
				sides[0].reserve((int)(ipts.size()*0.6));
				sides[1].reserve((int)(ipts.size()*0.6));

				int axis = 0;
				if (bbox.y_length() > bbox.x_length() && bbox.y_length() > bbox.z_length()) axis=1;
				if (bbox.z_length() > bbox.x_length() && bbox.z_length() > bbox.y_length()) axis=2;

				real_type split = bbox.centroid()[axis];

				for (unsigned i=0; i<ipts.size(); i++) {

					if (getpoint(ipts[i])[axis] < split)
						sides[0].push_back(ipts[i]);
					else
						sides[1].push_back(ipts[i]);
				}

				RecursiveTrim(sides[0], marked);
				RecursiveTrim(sides[1], marked);
			}
		}


		const GuidanceField &guidance;
		const GETPOINT &getpoint;
		const vector<real_type> &ideal;
    };


    protected:
    real_type rho;
    real_type min_step;
    real_type max_step;
    real_type reduction;

	private:
    real_type StepRequired(real_type dist, real_type ideal) const;


	// overlay stuff
    GetPointVector<Point3> o_kdGetPoint;
    vector<Point3> o_gfpoints;
    vector<real_type> o_ideal_lengths;
    typedef gtb::KDTree<int, real_type, GetPointVector<Point3> > o_kdtree_type;
    o_kdtree_type *o_kdtree;
	int o_lastkdsize;

    void* o_OrderedPointTraverseStart(const Point3 &p);
    bool o_OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void o_OrderedPointTraverseEnd(void *ctx);

	real_type o_MaxStepLength(const Point3 &p, real_type startlen);


	// we need to be able to add points to the overlay - can only do it when 
	// there aren't any traverses going on.
	volatile int o_nopt; // number of ordered point traverses currently going
	volatile bool o_pause;
	SLObject o_sl;
};




#endif

