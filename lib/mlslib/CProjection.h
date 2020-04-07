
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


#ifndef __CPROJECTION_2233_H
#define __CPROJECTION_2233_H

#include "mlslibdefs.h"
#include "knn.h"
#include "weight.h"
#include "poly.h"

#include <gtb/math/gradient.h>


MLSLIB_BEGIN_NAMESPACE

//
// Plane / polynomial fitting for a set of points
//

#define INFINITE_RADIUS ((REAL)1e8)

template <class REAL>
class CProjection
{
public:
    typedef REAL areal;

    typedef gtb::tsurfel_set<REAL> surfel_set;
    typedef gtb::tsurfelset_view<REAL> surfelset_view;
    typedef gtb::tPoint3<REAL> Point3;
    typedef gtb::tVector3<REAL> Vector3;
    typedef gtb::tPlane<REAL> Plane;
    typedef gtb::AMatc<REAL> AMatc;
    typedef gtb::AVec<REAL> AVec;
    typedef gtb::plane_transformation<REAL> plane_transformation;
    typedef gtb::ss_kdtree<surfel_set> kd_ss;

    typedef KNN<REAL> aKNN;
    typedef WeightFunction<REAL> WF;
    typedef Poly<REAL> lPoly;


    CProjection(
        surfel_set& points,
        int polynomial_degree);

    CProjection(
        surfel_set& points,
        int knn_radius,
        WF* theta,
        WF* radius_wf,
        int polynomial_degree,
        aKNN* knn);

    void print(); // print parameters

    void extract(const Point3& r, areal radius, surfelset_view& NN) const;
    void extract(const Point3& r, int K, surfelset_view& NN) const;

    //
    // the extract2 version are a 2-pass extract
    // first find the closest point x_r to r
    // and then query the neighborhood of x_r
    //
    void extract2(const Point3& r, areal radius, surfelset_view& NN) const;
    void extract2(const Point3& r, int K, surfelset_view& NN) const;

    //
    // The non-linear projection
    // Used to use powell....
    bool PowellProject(const Point3& r, Point3& r1, Vector3& n1) const;
    bool PowellProject(surfelset_view& nbhd, const Point3& r, Point3& r1, Vector3& n1, plane_transformation& T, lPoly* p, surfel_set& std_points) const;

    //
    // Approximated projection, i.e. TLS plane
    //
    bool approximated_projection(const Point3& r, Point3& r1, Vector3& n1);
    bool approximated_projection(surfelset_view& nbhd, const Point3& r, Point3& r1, Vector3& n1, plane_transformation& T, lPoly* p, surfel_set& std_points);


	//
	// the implicit function based projection from adamson/alexa's "approximating and intersecting surfaces from points"
	//

	// always use doubles for the gradients!
	// style parameter: 0=standard mls(always returns false!), 1=adamson with weighted ave of normals, 2=adamson with covariance normals from weighted average, 3=adamson with covariance normals from point

	typedef Gradient<double,3> Grad13;
	typedef Gradient<Grad13,3> Grad23;
	template <typename T>
	T adamson_eval(const surfelset_view &nbhd, const gtb::tPoint3<T> &p, REAL radius, int style) const;
    bool adamson_projection(const Point3& r, Point3& r1, Vector3& n1, int style, REAL *k1=NULL, REAL *k2=NULL) const;
    bool adamson_projection(surfelset_view& nbhd, const Point3& r, Point3& r1, Vector3& n1, int style, REAL *k1=NULL, REAL *k2=NULL) const;


    //
    // Shepard like interpolation
    // i.e. weighted average of nbhd
    //
    bool ShepardProject(const Point3& r, Point3& r1);
    bool ShepardProject(surfelset_view& nbhd, const Point3& r, Point3& r1, areal radius=INFINITE_RADIUS);
    
    //
    // Approximated plane and poly
    //
    static bool approximated_pnp(surfelset_view& nbhd, const Point3& r, WF* theta, plane_transformation& T, lPoly* p, surfel_set* std_points=0);
    static bool Powell_pnp(surfelset_view& nbhd, const Point3& r, WF* theta, plane_transformation& T, lPoly* p, surfel_set* std_points=0);

    static bool PlaneThroughPoint(surfelset_view& nbhd, const Point3& r, Plane& plane);
    static bool PlaneThroughCentroid(surfelset_view& nbhd, Plane& plane, Point3& centroid);
    static bool WeightedPlaneThroughPoint(surfelset_view& nbhd, const Point3& r, WF* theta, Plane& plane);

    static bool PowellMLSPlane(surfelset_view& nbhd, const Point3& r, areal radius, Plane& plane, WF* theta);


    static bool WeightedPolyFit(
        surfelset_view& nbhd, 
        plane_transformation& T, 
        const Point3& x, 
        WF* theta, 
        lPoly* poly,
        surfel_set& std_points);

    template <class CONT>
    static void PolyFit(CONT& std_points, lPoly* poly);
    static void PolyFit(surfelset_view& nbhd, plane_transformation& T, lPoly* poly, surfel_set& std_points);
    static void PolyFitWithResiduals(surfelset_view& nbhd, plane_transformation& T, lPoly* poly, AVec& residuals);

    static void cs2std(const surfelset_view& cs_nbhd, const plane_transformation& T, surfel_set& std_nbhd);
    static void std2cs(const surfelset_view& std_nbhd, const plane_transformation& T, surfel_set& cs_nbhd);

    //
    // T is either surfelset_view or surfelset_view
    //
    template <class T>
    static void compute_points_radius(gtb::ss_kdtree<T>& kd, int knn_radius);

    // Same, but using the KNN structure. (faster?)
    void compute_points_radius(aKNN* knn, unsigned knn_radius);

    // Initialize the radius_wf range
    void compute_radius_wf_factor();

    //
    // Compute the point radius for _points
    //
    void compute_points_radius();

    //
    // Compute the radius of a point
    // as the weighted average of points radiuses
    //
    areal point_radius(const Point3& x) const;
    areal point_radius(const surfelset_view& NN, const Point3& x) const;

    lPoly* gen_poly() const;
    static lPoly* gen_poly(int degree);

    int polynomial_degree() { return _polynomial_degree; }

    // compute the normal to the poly at 0, in world coordinates
    static void normal0(plane_transformation& T, lPoly* p, const Vector3& n, Vector3& n1);
    static void normal0(plane_transformation& T, lPoly* p, Vector3& n1);

    //
    void set_radius_factor(areal factor);

    kd_ss& get_kdtree();
    
    enum {MLS_MAX_ITERATIONS=200};

    surfel_set& get_points() { return _points; }
    const surfel_set& get_points() const { return _points; }
// protected:
    surfel_set& _points;
    kd_ss _kd;
    int _knn_radius; // # of KNN to extract when computing a points radius
    areal _radius_factor;
    aKNN* _knn; // A KNN structure: used by the extract methods


    WF* _theta;
    WF* _radius_wf;

    int _polynomial_degree;

protected:
    static bool PlaneInitialValues(surfelset_view& nbhd, const Point3& r, areal radius, areal& ir, areal& is, areal& it);

    //
    // A function object to evaluate the "value" of a plane.
    // used by the nonlinear optimization of the plane for the MLS projection.
    //
    //
    class PlaneValueObject
    {
    public:
        PlaneValueObject(
            const surfelset_view& nbhd,
            Plane& plane,
            const Point3& r,
            WF* theta
            );

        areal operator() (areal* t) ;
        areal derivative_value() const;
        void update_plane(Plane& plane);

        void get_q(Point3& q) const;
        const Point3& get_r() const;
        WF* get_theta();
        const surfelset_view& get_nbhd();
        Plane& get_plane();
    protected:
        const surfelset_view& _nbhd;
        const Point3& _r;
        Plane& _plane;
        WF* _theta;

		std::vector<areal> _weights;
    };


    //
    // A function object to evaluate the "value" of a plane defined using two angles and a distance
    // used by the nonlinear optimization of the plane for the MLS projection.
    //
    // Powell and the function it calls search from an initial point
    // by adding 1 (in each direction). But for our problem
    // 1 may be too large.
    // So... t is factored by some constant c*radius
    // And the angles are factored such that a change in angle
    // will not change "q" by more than the same factor.
    //
    class RSTPlaneValueObject
    {
    public:
        RSTPlaneValueObject(
            const surfelset_view& nbhd,
            Plane& plane,
            const Point3& r,
            WF* theta
            );
        areal operator() (areal* rst) const;
        void partial_derivatives(areal* rst, areal* pd) const;
        areal derivative_value() const { return _pvo.derivative_value(); }
        const PlaneValueObject& get_pvo() const { return _pvo; }

        mutable int evaluations;
        mutable int deriv_evaluations;
    protected:
        mutable PlaneValueObject _pvo;
    };
};

MLSLIB_END_NAMESPACE

#endif // __CPROJECTION_2233_H
