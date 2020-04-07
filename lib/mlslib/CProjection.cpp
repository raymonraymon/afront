
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


#include "stdafx.h"
#include <gtb/gtb.hpp>
#include <stat/stat.h>
#include "NR/nr.h"
#include "CProjection.h"
#include "poly.h"
#include "poly1.h"
#include "poly2.h"


//#include "PointsAnalysis.h"

//#define FIT_POLY_WITH_PLANE
//#include <viewer/debug_draw.h>

using gtb::aptr;
using gtb::ipow;
using gtb::nrran1;
using gtb::absv;
using gtb::infinit_bounded_decay;

MLSLIB_BEGIN_NAMESPACE

#define QUICK_NORMAL_ORIENTATION 1


template<class REAL>
CProjection<REAL>::CProjection(
    surfel_set& points,
    int polynomial_degree) :
        _points(points),
        _kd(points),
        _radius_factor(1.0),
        _knn(0),
        _polynomial_degree(polynomial_degree)
{
}

template<class REAL>
CProjection<REAL>::CProjection(
    surfel_set& points,
    int knn_radius,
    WF* theta,
    WF* radius_wf,
    int polynomial_degree,
    aKNN* knn
    ) :
        _points(points),
        _kd(points),
        _knn_radius(knn_radius),
        _radius_factor(1.0),
        _knn(knn),
        _theta(theta),
        _radius_wf(radius_wf),
        _polynomial_degree(polynomial_degree)

{
}

template<class REAL>
void CProjection<REAL>::extract(const Point3& r, areal radius, surfelset_view& NN) const
{
    _kd.extract(r, radius, NN);
}

template<class REAL>
void CProjection<REAL>::extract(const Point3& r, int K, surfelset_view& NN) const
{
    _kd.extract(r, (unsigned)K, NN);
}

template<class REAL>
void CProjection<REAL>::extract2(const Point3& r, areal radius, surfelset_view& NN) const
{
    int x_r_idx = _kd.tree->FindMin(r);
#if 0
    if (_knn)
    {
        _knn->extract(x_r_idx, radius, NN);
    }
    else
#endif
    {
        extract(_points.vertex(x_r_idx), radius, NN);
    }
}

template<class REAL>
void CProjection<REAL>::extract2(const Point3& r, int K, surfelset_view& NN) const
{
    int x_r_idx = _kd.tree->FindMin(r);
#if 0
    if (_knn)
    {
        _knn->extract(x_r_idx, (unsigned)K, NN);
    }
    else
#endif
    {
        extract(_points.vertex(x_r_idx), K, NN);
    }
}

template<class REAL>
bool CProjection<REAL>::PowellProject(const Point3& r, Point3& r1, Vector3& n1) const
{
    areal radius = point_radius(r); 

    surfelset_view nbhd(_points);
    extract2(r, radius, nbhd);
    plane_transformation T;
    aptr<lPoly> p(gen_poly());
    surfel_set std_points;

    return  PowellProject(nbhd, r, r1, n1, T, p, std_points);
}

/*
 * Input:
 *    nbhd - neighborhood,
 *    r - point to project
 *
 * Output:
 *    r1, n1 - the projected point and its normal
 *    T - the plane that was used
 *    p - the polynomial tha was used
 *    std_points - the neighborhood on the T coordinate space
 */
template<class REAL>
bool CProjection<REAL>::PowellProject(
    surfelset_view& nbhd, 
    const Point3& r, 
    Point3& r1, 
    Vector3& n1,
    plane_transformation& T,
    lPoly* p,
    surfel_set& std_points  
    ) const
{
    if (nbhd.size() < 6) return false;

    areal radius = point_radius(nbhd, r);
//    printf("radius: %f\n", radius);

//    aptr<WF> ltheta = _theta->clone();
    WF* ltheta = _theta;
    ltheta->set_radius(radius/3.0);

    Plane plane;
    if (!PowellMLSPlane(nbhd, r, radius, plane, ltheta))
    {
        return false;
    }

    T = plane_transformation(plane, r);

    Point3 q = T.FromPlane(Point3(0,0,0));

    WeightedPolyFit(nbhd, T, q, ltheta, p, std_points);

    Point3 r1_plane(0, 0, p->eval(0,0));
    r1 = T.FromPlane(r1_plane);

#if QUICK_NORMAL_ORIENTATION==1
    if ((nbhd.size() > 0) && nbhd.has_normals()) normal0(T, p, nbhd.normal(0), n1);
    else normal0(T, p, n1);
#else
    Vector3 dominant_n = dominant_orientation(nbhd);
    if (nbhd.has_normals()) normal0(T, p, dominant_n, n1);
#endif // QUICK_NORMAL_ORIENTATION

    return true;
}



// project using adamson/alexa's implicit function formulation
template<class REAL>
bool CProjection<REAL>::adamson_projection(const Point3& r, Point3& r1, Vector3& n1, int style, REAL *k1, REAL *k2) const {
    areal radius = point_radius(r); 

    surfelset_view nbhd(_points);
    extract2(r, radius, nbhd);
    return  adamson_projection(nbhd, r, r1, n1, style, k1, k2);
}


template<class REAL>
template <typename T>
T CProjection<REAL>::adamson_eval(const surfelset_view &nbhd, const gtb::tPoint3<T> &p, REAL radius, int style) const {

	REAL factor = -1.0 / (2*radius*radius);
	std::vector< gtb::tPoint3<T> > tpts(nbhd.size());	// just so we don't have to keep recreating this
	for (unsigned i=0; i<nbhd.size(); i++) {
		tpts[i] = gtb::tPoint3<T>((T)nbhd.vertex(i)[0], (T)nbhd.vertex(i)[1], (T)nbhd.vertex(i)[2]);
	}

	std::vector<T> weights(nbhd.size());
	T wsum=T(0);
	for (unsigned i=0; i<nbhd.size(); i++) {
		T distsq = gtb::tPoint3<T>::squared_distance(p, tpts[i]);
		weights[i] = exp(distsq * factor);
		wsum += weights[i];
	}

	T iwsum = 1 / wsum;

	// take weighted average of the points
	gtb::tPoint3<T> center((T)0,(T)0,(T)0);
	for (unsigned i=0; i<nbhd.size(); i++) {
		center.add_scaled(tpts[i], weights[i]);
	}
	center.scalar_scale(iwsum);


	// find the normal to use
	gtb::tVector3<T> norm((T)0,(T)0,(T)0);

	if (style > 1) {

		T cv[3] = { T(0), T(0), T(0) };

		for (unsigned i=0; i<nbhd.size(); i++) {
			gtb::tVector3<T> diff = tpts[i] - ((style==2) ? center : p);
			for (int j=0; j<3; j++)
				cv[j] += weights[i] * diff[j]*diff[j];
		}

//		for (int j=0; j<3; j++)
//			cv[j] = sqrt(cv[j]); // * iwsum);	// really need the iwsum scale??

		T **a = NR::amatrix<T>(1,3, 1,3);
		T **v = NR::amatrix<T>(1,3, 1,3);
		for (int r=0; r<3; r++) {
			for (int c=0; c<3; c++) {
				a[r+1][c+1] = cv[r]*cv[c];
				v[r+1][c+1] = 0;
			}
		}

		T d[3];
		int nrot;

		NR::jacobi<T>(a, 3, &d[-1], v, &nrot);

		int min = (fabs(d[0]) < fabs(d[1])) ? 0 : 1;
		if (fabs(d[2]) < fabs(d[min])) min = 2;

		norm[0] = v[0+1][min+1];
		norm[1] = v[1+1][min+1];
		norm[2] = v[2+1][min+1];

	} else {
		// weighted average of normals
		for (unsigned i=0; i<nbhd.size(); i++) {
			norm += weights[i] * gtb::tVector3<T>((T)nbhd.normal(i)[0], (T)nbhd.normal(i)[1], (T)nbhd.normal(i)[2]);
		}
		norm *= 1.0 / norm.length();	// since .normalize() checks for is_zero, which is undefined for gradient types
	}

	return norm.dot(p-center);
}


template<class REAL>
bool CProjection<REAL>::adamson_projection(surfelset_view& nbhd, const Point3& r, Point3& r1, Vector3& n1, int style, REAL *k1, REAL *k2) const {

	if (style==0 || (style==1 && !nbhd.has_normals())) return false;


	areal radius = point_radius(nbhd, r) / 2;

	r1 = r;

	int iter=0;
	while (1) {

		iter++;
		if (iter>50) {
			std::cerr<<"iteration failed to converge"<<std::endl;
			return false;
		}

		gtb::tPoint3<Grad13> p(r1[0], r1[1], r1[2]);
		p[0].gradient[0]=1;
		p[1].gradient[1]=1;
		p[2].gradient[2]=1;
		Grad13 e = adamson_eval(nbhd, p, radius, style);

		double step = -e.value / (e.gradient[0]*e.gradient[0] + e.gradient[1]*e.gradient[1] + e.gradient[2]*e.gradient[2]);

		areal vals[2];
		for (int i=0; i<2; i++) {

			int dir[2] = {1,-1};

			Point3 tr = r1;
			tr[0] += (areal)(dir[i]*step*e.gradient[0]);
			tr[1] += (areal)(dir[i]*step*e.gradient[1]);
			tr[2] += (areal)(dir[i]*step*e.gradient[2]);

			vals[i] = adamson_eval(nbhd, tr, radius, style);
		}

		// if (fabs(vals[0]) > fabs(vals[1]))
		// 	step = -step;

		r1[0] += (areal)(step*e.gradient[0]);
		r1[1] += (areal)(step*e.gradient[1]);
		r1[2] += (areal)(step*e.gradient[2]);

		n1[0]=(areal)e.gradient[0];
		n1[1]=(areal)e.gradient[1];
		n1[2]=(areal)e.gradient[2];

		// check for stopping
		if (step<1e-4 && fabs(e.value)<1e-1) {
			break;
		}

	}

	n1.normalize();


	// compute curvature if we want it
	if (k1 && k2) {
		gtb::tPoint3<Grad23> p((Grad13)r1[0], (Grad13)r1[1], (Grad13)r1[2]);
		p[0].gradient[0].value = p[0].value.gradient[0] = 1;
		p[1].gradient[1].value = p[1].value.gradient[1] = 1;
		p[2].gradient[2].value = p[2].value.gradient[2] = 1;
		Grad23 e = adamson_eval(nbhd, p, radius, style);


		typedef gtb::tmat3<double> mat3;

		mat3 H;
		H[0][0] = e.gradient[0].gradient[0];	H[0][1] = e.gradient[0].gradient[1];	H[0][2] = e.gradient[0].gradient[2];
		H[1][0] = e.gradient[1].gradient[0];	H[1][1] = e.gradient[1].gradient[1];	H[1][2] = e.gradient[1].gradient[2];
		H[2][0] = e.gradient[2].gradient[0];	H[2][1] = e.gradient[2].gradient[1];	H[2][2] = e.gradient[2].gradient[2];

		double gradient[3];
		gradient[0] = e.value.gradient[0];
		gradient[1] = e.value.gradient[1];
		gradient[2] = e.value.gradient[2];
	    
		gtb::tVector3<double> normal(gradient[0], gradient[1], gradient[2]);
		double gradmag = normal.length();
		normal /= gradmag;

		mat3 nnt(normal[0]*normal, normal[1]*normal, normal[2]*normal);
		mat3 I(gtb::tVector3<double>(1,0,0), gtb::tVector3<double>(0,1,0), gtb::tVector3<double>(0,0,1));
		mat3 P = I - nnt;

		mat3 G = (-1/gradmag) * P*H*P;


		double T = G[0][0] + G[1][1] + G[2][2];	// trace
		double F=0;	// frobenius norm
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				F += G[i][j]*G[i][j];
			}
		}
		F = sqrt(F);

		double st = sqrt(std::max(0.0, 2*F*F-T*T));

		*k1 = (areal)((T+st)/2);
		*k2 = (areal)((T-st)/2);

		if (gtb::absv(*k1) > gtb::absv(*k2)) std::swap(*k1, *k2);
	}

	return true;
}



template<class REAL>
bool CProjection<REAL>::approximated_projection(const Point3& r, Point3& r1, Vector3& n1)
{
    areal radius = point_radius(r);

    surfelset_view nbhd(_points);
    extract2(r,radius, nbhd);
    plane_transformation T;
    aptr<lPoly> p(gen_poly());
    surfel_set std_points;

    return  approximated_projection(nbhd, r, r1, n1, T, p, std_points);
}

template<class REAL>
bool CProjection<REAL>::approximated_projection (
    surfelset_view& nbhd, 
    const Point3& r, 
    Point3& r1, 
    Vector3& n1, 
    plane_transformation& T,
    lPoly* p,
    surfel_set& std_points
    )
{
    Point3 centroid = nbhd.centroid();

    areal radius = point_radius(nbhd, r);

//    aptr<WF> ltheta = _theta->clone();
    WF* ltheta = _theta;
    ltheta->set_radius(radius/3.0);

    Plane plane;
    if (!WeightedPlaneThroughPoint(nbhd, centroid, ltheta, plane))
    {
        return false;
    }
    T = plane_transformation(plane, centroid);

    WeightedPolyFit(nbhd, T, r, ltheta, p, std_points);

    Point3 local_r = T.ToPlane(r);
    Point3 r1_plane(local_r[0], local_r[1], p->eval(local_r[0], local_r[1]));
    r1 = T.FromPlane(r1_plane);

#if QUICK_NORMAL_ORIENTATION==1
    if (nbhd.has_normals()) normal0(T, p, nbhd.normal(0), n1);
#else
    Vector3 dominant_n = dominant_orientation(nbhd);
    if (nbhd.has_normals()) normal0(T, p, dominant_n, n1);
#endif // QUICK_NORMAL_ORIENTATION

    return true;
}

template<class REAL>
bool CProjection<REAL>::approximated_pnp(surfelset_view& nbhd, const Point3& r, WF* theta, plane_transformation& T, lPoly* p, surfel_set* std_points)
{
    Plane plane;
    if (!WeightedPlaneThroughPoint(nbhd, r, theta, plane)) return false;
    Point3 centroid = nbhd.centroid();
    T = plane_transformation(plane, centroid);

    if (std_points == 0)
    {
        surfel_set lstd_points;
        WeightedPolyFit(nbhd, T, r, theta, p, lstd_points);
    }
    else
    {
        WeightedPolyFit(nbhd, T, r, theta, p, *std_points);
    }

    return true;
}


template<class REAL>
bool CProjection<REAL>::ShepardProject(const Point3& r, Point3& r1)
{
    areal radius = point_radius(r);
    surfelset_view nbhd(_points);
    extract2(r,radius, nbhd);
    return ShepardProject(nbhd, r, r1, radius);
}

template<class REAL>
bool CProjection<REAL>::ShepardProject(surfelset_view& nbhd, const Point3& r, Point3& r1, areal radius)
{
    if (radius == INFINITE_RADIUS) 
    {
        radius = point_radius(r); 
    }
//    aptr<WF> ltheta = _theta->clone();
    WF* ltheta = _theta;
    ltheta->set_radius(radius/3.0);

    r1 = Point3::ZERO;
    areal weights = 0;
    int N = nbhd.size();
    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = nbhd.vertex(i);
        areal d2 = (r-xi).squared_length();
        areal wi = ltheta->weight(d2);
        r1.add_scaled(xi, wi);
        weights += wi;
    }
    r1.scalar_scale(((areal)1)/weights);
    return true;
}

template<class REAL>
bool CProjection<REAL>::Powell_pnp(surfelset_view& nbhd, const Point3& r, WF* theta, plane_transformation& T, lPoly* p, surfel_set* std_points)
{
    Plane plane;

    if (!PowellMLSPlane(nbhd, r, 0, plane, theta)) return false;

    T = plane_transformation(plane, r);

    if (std_points == 0)
    {
        surfel_set lstd_points;
        WeightedPolyFit(nbhd, T, r, theta, p, lstd_points);
    }
    else
    {
        WeightedPolyFit(nbhd, T, r, theta, p, *std_points);
    }

    return true;
}


template<class REAL>
bool CProjection<REAL>::PlaneThroughPoint(surfelset_view& nbhd, const Point3& r, Plane& plane)
{
    int n_points = nbhd.size();

    gtb::AMatc<REAL> MtM(3,3);
    MtM.set(0);

    for (int i = 0; i < n_points; ++i)
    {
        Vector3 pt(nbhd.vertex(i) - r);
        MtM(0,0) += pt[0] * pt[0];
        MtM(0,1) += pt[0] * pt[1];
        MtM(0,2) += pt[0] * pt[2];
        MtM(1,1) += pt[1] * pt[1];
        MtM(1,2) += pt[1] * pt[2];
        MtM(2,2) += pt[2] * pt[2];
    }

    AVec evalues(3);
    EigenUpperTriangular(MtM, evalues);


    //
    // I want the first eigen value to be much smaller than
    // other two.
    // this means
    // ev3 / ev2 < much_smaller_factor * ev2 / ev1
    // if all are positive (and they are)
    // Then ev3 * ev1 needs to be smaller than 2*ev2^2
    //
    if (evalues(0) * evalues(2) > 4*evalues(1)*evalues(1))
    {
        return false;
    }

    plane.normal() = Vector3(MtM(0,0), MtM(1,0), MtM(2,0));
    plane.setD(0);
    return true;
}

template<class REAL>
bool CProjection<REAL>::PlaneThroughCentroid(surfelset_view& nbhd, Plane& plane, Point3& centroid)
{
    centroid = nbhd.centroid();
    return PlaneThroughPoint(nbhd, centroid, plane);
}

template<class REAL>
bool CProjection<REAL>::WeightedPlaneThroughPoint(surfelset_view& nbhd, const Point3& r, WF* theta, Plane& plane)
{
    int n_points = nbhd.size();

    /*--------------------- max weight --------------------*/
    // first compute the maximal weight
    areal maxw = 0;
    for (int i = 0; i < n_points; ++i)
    {
        Vector3 pt(nbhd.vertex(i) - r);
        areal d2 = pt.squared_length();
        areal w = theta->weight(d2);
        if (maxw < w) maxw = w;
    }
    areal invmaxw = 1.0 / maxw;
    /*--------------------- [ max weight ] --------------------*/

    AMatc MtM(3,3);
    MtM.set(0);

    for (int i = 0; i < n_points; ++i)
    {
        Vector3 pt(nbhd.vertex(i) - r);
        areal d2 = pt.squared_length();
        areal w = theta->weight(d2) * invmaxw;

        MtM(0,0) += pt[0] * pt[0] * w;
        MtM(0,1) += pt[0] * pt[1] * w;
        MtM(0,2) += pt[0] * pt[2] * w;
        MtM(1,1) += pt[1] * pt[1] * w;
        MtM(1,2) += pt[1] * pt[2] * w;
        MtM(2,2) += pt[2] * pt[2] * w;
    }

    AVec evalues(3);
    EigenUpperTriangular(MtM, evalues);


    //
    // I want the first eigen value to be much smaller than
    // other two.
    // this means
    // ev3 / ev2 < much_smaller_factor * ev2 / ev1
    // if all are positive (and they are)
    // Then ev3 * ev1 needs to be smaller than 2*ev2^2
    //
    if (evalues(0) * evalues(2) > 4*evalues(1)*evalues(1))
    {
        return false;
    }

    plane.normal() = Vector3(MtM(0,0), MtM(1,0), MtM(2,0));
    plane.setD(0);
    return true;
}

#if 1
template<class REAL>
template<class CONT>
void CProjection<REAL>::PolyFit(CONT& std_points, lPoly* poly)
{
    int N = std_points.size();

    int NCF = poly->num_coefficientes();
    int DEG = poly->degree();

    AMatc A(NCF, NCF);
    AVec b(NCF);
    A.set(0);
    b.set(0);

    aptr<areal> cc(new areal[NCF]);

    /*--------- Prepare the matrix A and vector b ----------*/
    for (int k = 0; k < N; ++k)
    {
        const Point3& xi = std_points.vertex(k);


        //
        // Speedup: 
        //   - if our function is a*f_a(x,y) + b*f_b(x,y) + c*f_c(x,y) +...
        //     compute f_a(x,y), f_b(x,y), ...
        //     i.e. computes the "coefficients" of a,b,c...
        //   - Polynoms (old) pre-compute X^i, y^i, i=0:deg
        //
        {
            int cidx = 0;
            for (int p = 0; p <= DEG; ++p)
            {
                for (int ypwr = 0; ypwr <= p; ++ypwr)
                {
                    int xpwr = p - ypwr;
                    cc[cidx] = ipow(xi[0], xpwr) * ipow(xi[1], ypwr);
                    ++cidx;
                }
            }
        }

        areal *pr, *pc;
        for (int i = 0; i < NCF; ++i)
        {
            pc = pr = A.flatten() + i + i*NCF;
            for (int j = i; j < NCF; ++j)
            {
                *pr = *pc += cc[i]*cc[j]; // *w;
                pr += NCF;
                ++pc;
            }

            // xi[2] is the height of the function above the 
            // plane. (i.e. y_i)
            b[i] += xi[2] * cc[i]; // * w;
        }
    }

    /*--------- Solve Ax=b ----------*/
    {
        AVec x(NCF);
       
        SVDSolve(A, b, x);

        /*
         * Copy result to the coeeficients of the polynom
         */
        poly->SetCoefficients(x);
    }
}
#endif

#if 0
template void CProjection<float>::PolyFit(CProjection<float>::surfel_set&, Poly<float>*);
template void CProjection<double>::PolyFit(CProjection<double>::surfel_set&, Poly<double>*);

template void CProjection<float>::PolyFit(CProjection<float>::surfelset_view&, Poly<float>*);
template void CProjection<double>::PolyFit(CProjection<double>::surfelset_view&, Poly<double>*);
#endif

// std_points is output
template<class REAL>
void CProjection<REAL>::PolyFit(surfelset_view& nbhd, plane_transformation& T, lPoly* poly, surfel_set& std_points)
{
    cs2std(nbhd, T, std_points);
    PolyFit(std_points, poly);
}

//
// Fit a polynomial and return the residuals as well.
//
template<class REAL>
void CProjection<REAL>::PolyFitWithResiduals(surfelset_view& nbhd, plane_transformation& T, lPoly* poly, AVec& residuals)
{
    int N = nbhd.size();
    surfel_set std_points;
    PolyFit(nbhd, T, poly, std_points);

    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = std_points.vertex(i);
        areal h = poly->eval(xi[0], xi[1]);
        areal ri = xi[2] - h;
        residuals(i) = ri;
    }
}

template<class REAL>
bool CProjection<REAL>::WeightedPolyFit(
    surfelset_view& nbhd, 
    plane_transformation& T, 
    const Point3& x, 
    WF* theta, 
    lPoly* poly,
    surfel_set& std_points  // Output the std_points are stored here
                            // Some applications need that
    ) 
{
    int N = nbhd.size();
    cs2std(nbhd, T, std_points);

    int NCF = poly->num_coefficientes();
    int DEG = poly->degree();

    AMatc A(NCF, NCF);
    AVec b(NCF);
    A.set(0);
    b.set(0);

    aptr<areal> cc(new areal[NCF]);

    /*--------------------- max weight --------------------*/
    // first compute the maximal weight
    areal maxw = 0;
    areal minw = 1e9;
    for (int k = 0; k < N; ++k)
    {
        // unused const Point3& xi = std_points.vertex(k);

        areal d = (nbhd.vertex(k) - x).squared_length();
        areal w = theta->weight(d);
        if (maxw < w) maxw = w;
        if (minw > w) minw = w;
    }
    if (maxw == 0)
    {
//        DebugBreak();
        return false;
    }
    areal invmaxw = 1.0 / maxw;
    /*--------------------- [ max weight ] --------------------*/

    /*--------- Prepare the matrix A and vector b ----------*/
    for (int k = 0; k < N; ++k)
    {
        const Point3& xi = std_points.vertex(k);

        areal d = (nbhd.vertex(k) - x).squared_length();
        areal w = theta->weight(d) * invmaxw;

        //
        // Speedup?
        //   - if our function is a*f_a(x,y) + b*f_b(x,y) + c*f_c(x,y) +...
        //     compute f_a(x,y), f_b(x,y), ...
        //     i.e. computes the "coefficients" of a,b,c...
        //   - Polynoms (old) pre-compute X^i, y^i, i=0:deg
        //
        {
            int cidx = 0;
            for (int p = 0; p <= DEG; ++p)
            {
                for (int ypwr = 0; ypwr <= p; ++ypwr)
                {
                    int xpwr = p - ypwr;
                    cc[cidx] = ipow(xi[0], xpwr) * ipow(xi[1], ypwr);
                    ++cidx;
                }
            }
        }

        areal *pr, *pc;
        for (int i = 0; i < NCF; ++i)
        {
            pc = pr = A.flatten() + i + i*NCF;
            for (int j = i; j < NCF; ++j)
            {
                *pr = *pc += cc[i]*cc[j] *w;
                pr += NCF;
                ++pc;
            }

            // xi[2] is the height of the function above the 
            // plane. (i.e. y_i)
            b[i] += xi[2] * cc[i] * w;
        }
    }

    /*--------- Solve Ax=b ----------*/
    {
        if (!isfinite(A(0,0)))
        {
//            DebugBreak();
            return false;
        }
        else
        {
            AVec x(NCF);
           
            SVDSolve(A, b, x);

            /*
             * Copy result to the coeeficients of the polynom
             */
            poly->SetCoefficients(x);
        }
    }
    return true;
}

template<class REAL>
void CProjection<REAL>::cs2std(const surfelset_view& world_nbhd, const plane_transformation& T, surfel_set& std_nbhd)
{
    int N = world_nbhd.size();
    std_nbhd.clear();
    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = world_nbhd.vertex(i);
        std_nbhd.insert_vertex(T.ToPlane(xi));
    }
}

template<class REAL>
void CProjection<REAL>::std2cs(const surfelset_view& std_nbhd, const plane_transformation& T, surfel_set& world_nbhd)
{
    int N = std_nbhd.size();
    world_nbhd.clear();
    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = std_nbhd.vertex(i);
        world_nbhd.insert_vertex(T.FromPlane(xi));
    }
}

template<class REAL>
void CProjection<REAL>::compute_points_radius()
{
    if (!_points.has_radius())
    {
        if (_knn) compute_points_radius(_knn, _knn_radius);
        else compute_points_radius(_kd, _knn_radius);
    }

    compute_radius_wf_factor();
}

template<class REAL>
void CProjection<REAL>::compute_radius_wf_factor()
{
    // Set the radius weight function parameter
    // to the median of radiuses.
    std::vector<areal> sample_radiuses;
    static const int NSAMPLES=100;
    sample_radiuses.reserve(NSAMPLES);
    for (int i = 0; i < NSAMPLES; ++i)
    {
        sample_radiuses.push_back(_points.radius(nrran1()%_points.size()));
    }
    areal medr;
    stattools::median(sample_radiuses.begin(), sample_radiuses.end(), medr);
    _radius_wf->set_radius(medr);
}

//
// Compute the radius of each point
// that is the distance of the Kth nearest neighbor
//
template<class REAL>
template <class T>
void CProjection<REAL>::compute_points_radius(gtb::ss_kdtree<T>& kd, int knn_radius)
{
    T& points = kd.get_points();
    points.clear_radius();
    int N = points.size();

    for (int i = 0; i < N; ++i)
    {
        const Point3& x = points.vertex(i);
        surfelset_view nbhd(points);
        kd.extract(x, (unsigned)knn_radius, nbhd);
        areal radius = (x-nbhd.vertex(knn_radius-1)).length();
        points.insert_radius(radius);
    }
}

template<class REAL>
void CProjection<REAL>::compute_points_radius(aKNN* knn, unsigned knn_radius)
{
    surfel_set& points = knn->get_points();

    points.clear_radius();
    unsigned N = points.size();
    surfelset_view nbhd(points);
    for (unsigned i = 0; i < N; ++i)
    {
        nbhd.clear();
        const Point3& xi = points.vertex(i);
        knn->extract(i, knn_radius, nbhd);
        unsigned Ns = nbhd.size();
        if (Ns < knn_radius)
        {
            printf("wrong nbhd size (%d) for %d %g %g %g\n", Ns, i, xi[0], xi[1], xi[2]);
        }
        areal radius = (xi-nbhd.vertex(Ns-1)).length();
        points.insert_radius(radius);
    }
}

template<class REAL>
typename CProjection<REAL>::areal CProjection<REAL>::point_radius(const Point3& x) const
{
    surfelset_view NN(_points);
    extract2(x, _knn_radius, NN);
    return point_radius(NN, x);
}

template<class REAL>
typename CProjection<REAL>::areal CProjection<REAL>::point_radius(const surfelset_view& NN, const Point3& x) const
{
    areal radius = 0;
    areal normalizer = 0;
    int n_nbrs = NN.size();
    for (int i = 0; i < n_nbrs; ++i)
    {
        areal d2 = (x-NN.vertex(i)).squared_length();
        areal w = _radius_wf->weight(d2);
        radius += NN.radius(i)*w;
        normalizer += w;
    }
    if (absv(normalizer)<1e-8)
    {
        printf("r:%f (%d)\n", radius, NN.size());
        printf("Point with bad radius\n");
        return 0;
        radius = NN.radius(0);
    }
    else radius /= normalizer;
    return radius * _radius_factor;
}

template<class REAL>
typename CProjection<REAL>::lPoly* CProjection<REAL>::gen_poly(int degree)
{
    switch(degree)
    {
    case 1:
        return new Poly1<REAL>;
    case 2:
        return new Poly2<REAL>;
    default:
        assert(0);
        return 0; // warning, not all control paths...
    }
}

template<class REAL>
typename CProjection<REAL>::lPoly* CProjection<REAL>::gen_poly() const
{
    return gen_poly(_polynomial_degree);
}

template<class REAL>
void CProjection<REAL>::normal0(plane_transformation& T, lPoly* p, const Vector3& n, Vector3& n1)
{
    p->Normal0(n1);
    n1 = T.PlaneNormal2World(n1);
    if (n1.dot(n) < 0) n1.flip();
}

template<class REAL>
void CProjection<REAL>::normal0(plane_transformation& T, lPoly* p, Vector3& n1)
{
    p->Normal0(n1);
    n1 = T.PlaneNormal2World(n1);
}

/*------------------ PlaneValueObject ------------------*/
template<class REAL>
CProjection<REAL>::PlaneValueObject::PlaneValueObject(
      const surfelset_view& nbhd,
      Plane& plane,
      const Point3& r,
      WF* theta
      ) :
    _nbhd(nbhd),
    _r(r),
    _plane(plane),
    _theta(theta),
	_weights(nbhd.size())
{
}

template<class REAL>
typename CProjection<REAL>::areal CProjection<REAL>::PlaneValueObject::operator() (areal* t) 
{
    const Vector3& n = _plane.normal();

    areal normalizer = 0.0;
    areal value = 0.0;

    Point3 q = _r + n*t[1];

    unsigned N = _nbhd.size();
#if 0
    areal sumabsxiwi=0;

    // Faster, but numerically unstable
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& x = _nbhd.vertex(i);
        Vector3 qtox = x-q;
        areal qdist = n.dot(qtox);
        areal qdist2 = qdist * qdist;
        areal weight = _theta->weight(qtox.squared_length());

        normalizer += weight;
        value += qdist2 * weight;
        sumabsxiwi += gtb::absv(qdist2) * weight;
    }
//    assert(normalizer > 0);
    value /= normalizer;
    printf("Cond: %f\n", sumabsxiwi / value);
//    value /= (normalizer*normalizer);
#else // 0/1


	// What is faster? recompute in 2 loops or save the result?
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& x = _nbhd.vertex(i);
        Vector3 qtox = x-q;
        areal weight = _theta->weight(qtox.squared_length());
		_weights[i] = weight;
        normalizer += weight;
    }

    if (normalizer < 1e-8) return 1e8;
    normalizer = 1.0 / normalizer;




#ifdef FIT_POLY_WITH_PLANE


	int min = fabs(n[0])<fabs(n[1]) ? 0 : 1;
	min = fabs(n[2])<fabs(n[min]) ? 2 : min;
    Vector3 udir(0,0,0);
	udir[min] = 1;
	udir = n.cross(udir);
	udir.normalize();
	Vector3 vdir = udir.cross(n);


	gtb::AMat<double> A(N, 5);
	gtb::AVec<double> b(N);
	gtb::AVec<double> x(5);

	for (unsigned i=0; i<N; ++i) {

        const Point3& x = _nbhd.vertex(i);
        Vector3 qtox = x-q;
		
		areal u = udir.dot(qtox);
		areal v = vdir.dot(qtox);

        areal weight = _theta->weight(qtox.squared_length()) * normalizer;

		A(i,0) = u*u * weight;
		A(i,1) = u*v * weight;
		A(i,2) = v*v * weight;
		A(i,3) = u * weight;
		A(i,4) = v * weight;
//		A(i,5) = 1 * weight;
		b(i) = n.dot(qtox) * weight;
	}

	gtb::SVDSolve(A, b, x);

//	std::cerr<<x<<std::endl;


	for (unsigned i=0; i<N; ++i) {

        const Point3& px = _nbhd.vertex(i);
        Vector3 qtox = px-q;
		
		areal u = udir.dot(qtox);
		areal v = vdir.dot(qtox);


		areal polyheight = x[0]*u*u + x[1]*u*v + x[2]*v*v + x[3]*u + x[4]*v/* + x[5]*/;
		areal pointheight = n.dot(qtox);
		areal dist = pointheight - polyheight;

        areal weight = _theta->weight(qtox.squared_length()) * normalizer;
		value += dist*dist*weight;

/*
		std::cerr<<"u: "<<u<<std::endl;
		std::cerr<<"v: "<<v<<std::endl;
		std::cerr<<"pointheight: "<<pointheight<<std::endl;
		std::cerr<<"polyheight: "<<polyheight<<std::endl;
		std::cerr<<"weight: "<<weight<<std::endl;
		std::cerr<<"dist: "<<dist<<std::endl;
		std::cerr<<"wdist: "<<(dist*dist*weight)<<std::endl;
		std::cerr<<std::endl;
*/
	}


	if (0){
		const int num=25;
		const areal step = (areal)0.1;


		dbgClear();


		// the plane
		{
			std::vector<Point3> plane(4);
			plane[0] = q + -num*step*udir + -num*step*vdir;
			plane[1] = q +  num*step*udir + -num*step*vdir;
			plane[2] = q +  num*step*udir +  num*step*vdir;
			plane[3] = q + -num*step*udir +  num*step*vdir;

			std::vector<Point3> pcolor(4);
			pcolor[0] = Point3(0,0,1);
			pcolor[1] = Point3(0,0,1);
			pcolor[2] = Point3(0,0,1);
			pcolor[3] = Point3(0,0,1);

			DbgPoly::add(0, plane, pcolor);

		}


		// the coordinate system
		{
			std::vector<Point3> line(2);
			line[0] = q;

			line[1] = q + num*step*udir;
			DbgPLines::add(0, line, 3);

			line[1] = q - num*step*udir;
			DbgPLines::add(0, line, 3);

			line[1] = q + num*step*vdir;
			DbgPLines::add(0, line, 3);

			line[1] = q - num*step*vdir;
			DbgPLines::add(0, line, 3);
		}



		std::vector<Point3> tri(3);

		for (int yi=-num; yi<num; yi++) {

			areal ry = yi*step;

			for (int xi=-num; xi<num; xi++) {

				areal rx = xi*step;


				DbgPoly::add(0, q+rx*udir+ry*vdir + ((areal)(x[0]*rx*rx + x[1]*rx*ry + x[2]*ry*ry/* + x[3]*rx + x[4]*ry + x[5]*/))*n, Point3(0,1,0),
								q+(rx+step)*udir+ry*vdir + ((areal)(x[0]*(rx+step)*(rx+step) + x[1]*(rx+step)*ry + x[2]*ry*ry/* + x[3]*(rx+step) + x[4]*ry + x[5]*/))*n, Point3(0,1,0),
								q+rx*udir+(ry+step)*vdir + ((areal)(x[0]*rx*rx + x[1]*rx*(ry+step) + x[2]*(ry+step)*(ry+step)/* + x[3]*rx + x[4]*(ry+step) + x[5]*/))*n, Point3(0,1,0));

				DbgPoly::add(0, q+(rx+step)*udir+ry*vdir + ((areal)(x[0]*(rx+step)*(rx+step) + x[1]*(rx+step)*ry + x[2]*ry*ry/* + x[3]*(rx+step) + x[4]*ry + x[5]*/))*n, Point3(0,1,0),
								q+rx*udir+(ry+step)*vdir + ((areal)(x[0]*rx*rx + x[1]*rx*(ry+step) + x[2]*(ry+step)*(ry+step)/* + x[3]*rx + x[4]*(ry+step) + x[5]*/))*n, Point3(0,1,0),
								q+(rx+step)*udir+(ry+step)*vdir + ((areal)(x[0]*(rx+step)*(rx+step) + x[1]*(rx+step)*(ry+step) + x[2]*(ry+step)*(ry+step)/* + x[3]*(rx+step) + x[4]*(ry+step) + x[5]*/))*n, Point3(0,1,0));


			}
		}


		std::cerr<<"value:"<<value<<std::endl;
		extern void redrawAndWait(int key, bool force);
		redrawAndWait(0,false);
	}



#else

	areal qdot = n.dot(q);
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& x = _nbhd.vertex(i);
//        Vector3 qtox = x-q;
        areal qdist = n.dot(x)-qdot; //n.dot(qtox);
        areal qdist2 = qdist * qdist;
        areal weight = _weights[i] * normalizer; //_theta->weight(qtox.squared_length()) * normalizer;

        value += qdist2 * weight;
    }
#endif

#endif
//    printf("%0.20f\n", value);
    return value;
}

/*
 * compute the derivative of the function with respect to t at the current t
 * for now, derive the function without the normalization factor
 * Assuming gaussian weight function
 */
template<class REAL>
typename CProjection<REAL>::areal CProjection<REAL>::PlaneValueObject::derivative_value() const
{
    const Vector3& n = _plane.normal();
    areal t = _plane.D();

    areal value = 0.0;
    Point3 q = _r + n*t;
    GaussianWeightFunction<REAL>* p_gt = static_cast<GaussianWeightFunction<REAL>*>(_theta); // HACK!
//    areal sigma2 = 1.0 / (p_gt->_sigma * p_gt->_sigma);
    areal sigma2 = 1.0 / p_gt->_sigma;

    unsigned N = _nbhd.size();
    for (unsigned i = 0; i < N; ++i)
    {
        const Point3& x = _nbhd.vertex(i);
        Vector3 qtox = x-q;
//        Vector3 qtor = x-_r;
        areal qdist = n.dot(qtox);
        areal qdist2 = qdist * qdist;
//        areal rdist = n.dot(qtor);
        areal weight = _theta->weight(qtox.squared_length());

//        value += (t - rdist)* (1.0 + qdist2) * weight;
        value += qdist * (1.0 - qdist2 * sigma2) * weight;
    }
    return 2.0 * value;
}

template<class REAL>
void CProjection<REAL>::PlaneValueObject::update_plane(Plane& plane)
{
    _plane = plane;
}

template<class REAL>
void CProjection<REAL>::PlaneValueObject::get_q(Point3& q) const
{
    q = _r + _plane.normal() * _plane.D();
}

template<class REAL>
const typename CProjection<REAL>::Point3& CProjection<REAL>::PlaneValueObject::get_r() const
{
    return _r;
}

template<class REAL>
typename CProjection<REAL>::WF* CProjection<REAL>::PlaneValueObject::get_theta()
{
    return _theta;
}

template<class REAL>
const typename CProjection<REAL>::surfelset_view& CProjection<REAL>::PlaneValueObject::get_nbhd()
{
    return _nbhd;
}

template<class REAL>
typename CProjection<REAL>::Plane& CProjection<REAL>::PlaneValueObject::get_plane()
{
    return _plane;
}

/*------------------ RSTPlaneValueObject ------------------*/

template<class REAL>
CProjection<REAL>::RSTPlaneValueObject::RSTPlaneValueObject(
            const surfelset_view& nbhd,
            Plane& plane,
            const Point3& r,
            WF* theta
            ) : 
    evaluations(0),
    deriv_evaluations(0),
    _pvo(nbhd, plane, r, theta)
{
}


template<class REAL>
typename CProjection<REAL>::areal CProjection<REAL>::RSTPlaneValueObject::operator() (areal* rst) const
{
    static areal t0[2] = {0.0, 0.0};

    {
	gtb::tCPlaneRST<REAL> p(rst[1], rst[2], rst[3]);
	Plane temp(p);
    _pvo.update_plane(temp);
    }
    t0[1] = rst[3];
    areal v = _pvo(t0);
   
//    ++evaluations;
    return v;
}

template<class REAL>
void CProjection<REAL>::RSTPlaneValueObject::partial_derivatives(areal* rst, areal* pd) const
{
    static areal t0[2] = {0.0, 0.0};

    t0[1] = rst[3];
    {
	gtb::tCPlaneRST<REAL> p(rst[1], rst[2], t0[1]);
	Plane temp(p);
    _pvo.update_plane(temp);
    }

#if 1
    // the value of the function at the input point
    areal vp = _pvo(t0);

    //
    // compute partial derivatives numerically with h step
    // that is (f(p+h)-f(p))/h
    //
//    areal h = _tolerance*0.1;
//    areal h = max3(gtb::absv(rst[1]), gtb::absv(rst[2]),gtb::absv(rst[3]))*0.01;
//    printf("%f\n", h);
    areal h1,h2,h3;
#if 0
    h1=h2=max2(absv(rst[1]), absv(rst[2]))*(areal)0.1;
    h3=absv(rst[3])*(areal)0.1;
#else
    h1 = h2 = h3 = gtb::fp_sqrt_precision<REAL>();
#endif

    t0[1] = rst[3]+h3;
    areal vp_th = _pvo(t0);
    t0[1] = rst[3];
    {
	gtb::tCPlaneRST<REAL> p(rst[1]+h1, rst[2], t0[1]);
	Plane temp(p);
    _pvo.update_plane(temp);
    }
    areal vp_rh = _pvo(t0);
    {
	gtb::tCPlaneRST<REAL> p(rst[1], rst[2]+h1, t0[1]);
	Plane temp(p);
    _pvo.update_plane(temp);
    }
    areal vp_sh = _pvo(t0);

    pd[1] = (vp_rh - vp)/h1;
    pd[2] = (vp_sh - vp)/h2;
    pd[3] = (vp_th - vp)/h3;
#else
    //
    // A Higher accuracy numerical derivative:
    // df = (f(x+h) - f(x-h))/ (2h)
    // Numerical mathematics and computing / Cehney, pg 100
    areal h1,h2,h3;
    h1=h2=max2(absv(rst[1]), absv(rst[2]))*(areal)0.1;
    h3=absv(rst[3])*(areal)0.1;

    t0[1] = rst[3]+h3;
    areal vpt_ph = _pvo(t0);
    t0[1] = rst[3]-h3;
    areal vpt_mh = _pvo(t0);
    t0[1] = rst[3];

    {
	Plane temp(CPlaneRST(rst[1]+h1, rst[2], t0[1]));
    _pvo.update_plane(temp);
    }
    areal vpr_ph = _pvo(t0);

    {
	Plane temp(CPlaneRST(rst[1]-h1, rst[2], t0[1]));
    _pvo.update_plane(temp);
    }
    areal vpr_mh = _pvo(t0);

    {
	Plane temp(CPlaneRST(rst[1], rst[2]+h1, t0[1]));
    _pvo.update_plane(temp);
    areal vps_ph = _pvo(t0);

    {
	Plane temp(CPlaneRST(rst[1], rst[2]-h1, t0[1]));
    _pvo.update_plane(temp);
    }
    areal vps_mh = _pvo(t0);

    pd[1] = (vpr_ph-vpr_mh)/(2*h1);
    pd[2] = (vps_ph-vps_mh)/(2*h1);
    pd[3] = (vpt_ph-vpt_mh)/(2*h3);
#endif

#if 0
    printf("PD: at <%g %g %g>: <%g %g %g> (%g %g %g)\n", 
        rad_to_deg(rst[1]), rad_to_deg(rst[2]), rst[3],
        pd[1], pd[2], pd[3],
        h1,h2,h3);
#endif
//    ++deriv_evaluations;
}

/*
 * Compute the initial values for r,s,t for the MLS plane
 * Used by the optimization process
 *
 * Input:
 *    nbhd
 *    r - projected point
 *    radius - of neighborhood
 * Ourput:
 *    ir, is, it - initial values for r,s,t
 *
 * Return false if failed
 */
template<class REAL>
bool CProjection<REAL>::PlaneInitialValues(surfelset_view& nbhd, const Point3& r, areal radius, areal& ir, areal& is, areal& it)
{
    Plane plane;
    Point3 center;
    if (!PlaneThroughCentroid(nbhd, plane ,center)) return false;

    //
    // First estimation:
    //    project r on the plane
    //
    Vector3 n1 = plane.normal();

    // Project r on the plane
    Vector3 r_to_center = r - center;
    areal it1 = r_to_center.dot(n1);
    if (it1 > 0)
    {
        n1.flip();
    }
    else
    {
        it1 = -it1;
    }

    //
    // second estimation:
    //   q=center
    //
    Vector3 n2 = r_to_center;
    areal it2 = r_to_center.length();
    n2 /= it2;
    it2 = -it2;
    if (n2.dot(n1))
    {
        n2.flip();
        it2 = -it2;
    }

    //
    // Compute the weights of the two estimations
    //
    areal l = gtb::absv(it2);
    areal alpha = infinit_bounded_decay(l, radius, radius*(areal)0.1);
    Vector3 n = n1*alpha + n2*(((areal)1) - alpha);
    n.normalize();
    it = it1 * alpha + it2 * (((areal)1) - alpha);

#if 0
    printf(
        "Initials: <%g (%g %g %g)>*%g + <%g (%g %g %g)>*%g\n",
        it1, n1[0], n1[1], n1[2], alpha,
        it2, n2[0], n2[1], n2[2], 1.0 - alpha);
#endif



//    Point3 q = r - n * it;

    // compute ir, is
    is = asin(n[2]);
    ir = atan2(n[1], n[0]);
    // gtb::euclidian_to_spherical(n, ir, is); // This just won't get instantiated :(

        /*
         * Verified in debug session that the ppp.normal === n
            CPlaneRST rstptest(initial_r, initial_s, initial_t);
            Plane ppp(rstptest);
            areal xxx = 2;
         */

    
    /* // DEBUG the initial values
    plane = CPlaneRST(ir, is, it); 
    if(1 - plane.normal().dot(n) > 1e-5) printf("Bad plane\n");
    / **/
    return true;
}

//
// Compute the optimized plane to a point for
// the MLS projection procedure
//
#define OLDPOWELL 0
template<class REAL>
bool CProjection<REAL>::PowellMLSPlane(surfelset_view& nbhd, const Point3& r, areal radius, Plane& plane, WF* theta)
{
    areal initial_r, initial_s, initial_t;
    PlaneInitialValues(nbhd, r, theta->influence_radius(0.95), initial_r, initial_s, initial_t);

    RSTPlaneValueObject pvo(nbhd, plane, r, theta);

    areal* rst = NR::avector<areal>(1,3); // same as R,S,T in CPlaneRST
    rst[1] = initial_r;
    rst[2] = initial_s;
    rst[3] = initial_t;

//    energy_visualizer("e_mlsplane.obj", pvo, rst);

#if OLDPOWELL==1
    areal **xi = NR::amatrix<areal>(1, 3, 1, 3);
    xi[1][1] = 1; xi[1][2] = 0; xi[1][3] = 0;
    xi[2][1] = 0; xi[2][2] = 1; xi[2][3] = 0;
    xi[3][1] = 0; xi[3][2] = 0; xi[3][3] = 1;
#endif

    areal ftol = 1e-3f;        // BUGBUG should be a parameter
    int iter;                   // Placeholder for the #iterations taken.
    areal fret;                // Function value at "converged" point

//    plane = Plane(CPlaneRST(initial_r, initial_s, initial_t)); Debugging initial value

#if OLDPOWELL==1
    if (!NR::powell(rst, xi, 3, ftol, &iter, &fret, pvo))
#else
    const_mem_fun2_t<void, RSTPlaneValueObject, areal*, areal*> pvod(&pvo, &RSTPlaneValueObject::partial_derivatives);
    if (!NR::frprmn(rst, 3, ftol, &iter, &fret, pvo, pvod))
//    if (!NR::frprmn(rst, 3, ftol, &iter, &fret, pvo, mem_fun2<void,RSTPlaneValueObject>(&pvo, &RSTPlaneValueObject::partial_derivatives)))
//    if (!mygd(rst, 3, ftol, &iter, &fret, pvo, pvod))
#endif
    {
//        printf("CProjection<REAL>::PowellMLSPlane did not converge\n");
        return false;
    }
//    printf("converged in %d iterations, %d %d evaluations\n", iter, pvo.evaluations, pvo.deriv_evaluations);

#if 0
    {
        const PlaneValueObject& xpvo = pvo.get_pvo();
        areal t1[2] = {0, plane.D()-1e-6};
        areal t2[2] = {0, plane.D()+1e-6};
        areal numderiv = (xpvo(t2) - xpvo(t1))/(t2[1]- t1[1]);
        printf("dt = %g, numerically: %g\n", pvo.derivative_value(), numderiv);
    }
#endif

#if OLDPOWELL==1
    NR::free_amatrix(xi, 1, 3, 1, 3);
#endif

    NR::free_avector(rst, 1, 3);
    return true;
}


template<class REAL>
void CProjection<REAL>::set_radius_factor(areal factor)
{
    _radius_factor = factor;
}


template<class REAL>
typename CProjection<REAL>::kd_ss& CProjection<REAL>::get_kdtree()
{
    return _kd;
}

template<class REAL>
void CProjection<REAL>::print()
{
    printf(
        "CProjection:\n"
        "    Knn radius: %d\n"
        "    Radius factor: %f\n",
        _knn_radius,
        _radius_factor);
}

template class CProjection<float>;
template class CProjection<double>;

MLSLIB_END_NAMESPACE

