
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


#include <gtb/gtb.hpp>
#ifndef WIN32
#include <gtb/graphics/plane.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/math/mathtools.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/plane.ipp>
#undef inline
#endif

GTB_BEGIN_NAMESPACE



template<class T>
bool tPlane<T>::intersects(const tSegment3<T> &segment,
		       tPoint3<T> &hit_point,
		       value_type &hit_time) const
{
	// n (p + t.d) + D = 0
	// t = -(D + n.p) / n.d

	const tPoint3<T> &p = segment.source();
	tVector3<T> d = segment.direction();
	value_type nd = _n * d;
	value_type np = _n * p;

	if (nd == 0.0) {
		if ((np + _d) == 0.0) {
			// segment in the plane
			hit_point = p;
			hit_time = 0.0;
			return true;
		} else {
			// segment parallel to the plane
			return false;
		}
	} else {
		hit_time = -(_d + np) / nd;
		if ((hit_time >= 0.0) && (hit_time <= 1.0)) {
			hit_point = p + (hit_time * d);
			return true;
		} else {
			return false;
		}
	}
}

template<class T>
bool tPlane<T>::intersects(const tRay3<T> &ray,
			   tPoint3<T> &hit_point,
			   value_type &hit_time) const
{
	// n (p + t.d) + D = 0
	// t = -(D + n.p) / n.d

	const tPoint3<T> &o = ray.source();
	const tVector3<T> &d = ray.direction();
	value_type no = _n.dot(o);
	value_type nd = _n.dot(d);

	if (no == 0.0) {
		if ((no + _d) == 0.0) {
			// segment in the plane
			hit_point = o;
			hit_time = 0.0;
			return true;
		} else {
			// segment parallel to the plane
			return false;
		}
	} else {
		hit_time = -(no + _d) / nd;
		if (hit_time >= 0.0) {
			hit_point = o + (hit_time * d);
			return true;
		} else {
			return false;
		}
	}
}


//
// BUGBUG: assume value_type is double
//
template<class T>
void tPlane<T>::read(FILE *fp)
{
    double x;
    read_double(&x, fp); _n[0] = x;
    read_double(&x, fp); _n[1] = x;
    read_double(&x, fp); _n[2] = x;
    read_double(&x, fp); _d = x;
}

template<class T>
void tPlane<T>::write(FILE *fp) const
{
    write_double(_n[0], fp);
    write_double(_n[1], fp);
    write_double(_n[2], fp);
    write_double(_d, fp);
}


template<class T>
void tPlane<T>::flip_orientation()
{
    _n.flip();
    _d = -_d;
}

/*
 * Intesection of two planes forms a line
 * Assume non-parallel planes
 * Assume n is normalized
 */
template<class T>
bool tPlane<T>::intersect(const tPlane& P, tLine3<T>& l)
{
    const tVector3<T>& n1 = _n;
    const tVector3<T>& n2 = P._n;

    tVector3<T>& nl = l.direction();
    nl = n1.cross(n2);
    nl.normalize();

    // theoretically we would set l.origin.z = 0, but this may lead
    // to degenerecies, thus we set l.origin[idx3] = 0
    // and the other two indices are the remaining dimensions
    int idx1,idx2,idx3;
    max3(fabs(nl[0]), fabs(nl[1]), fabs(nl[2]), idx3);
    idx1 = (idx3+1)%3;
    idx2 = (idx1+1)%3;

    // solve the 2x2 matrix Ax=b, A=(a11,a12,a21,a22), b=(b1,b2)
    // A=[n1;n2]
    // b = [n1*p1; n2*p2]
    double a11,a12,a21,a22;
    double b1,b2;
    a11 = n1[idx1];
    a12 = n1[idx2];
    a21 = n2[idx1];
    a22 = n2[idx2];
    b1 = n1.dot(origin());
    b2 = n2.dot(P.origin());
    double det = a11*a22-a12*a21;

    if (fabs(det) < 1e-12) return false;

    tPoint3<T>& o = l.origin();
    o[idx1] = (a22*b1 - a12*b2)/det;
    o[idx2] = (a11*b2 - a21*b1)/det;
    o[idx3] = 0;

//    printf("%g %g\n", n1.dot(o-origin()), n2.dot(o-P.origin()));

    return true;
}

template<class T>
tPoint3<T> tPlane<T>::project(const tPoint3<T>& p) const
{
    tPoint3<T> o = origin();
    tPoint3<T> proj_p = p - _n.dot(p-o)*_n;

    assert(distance(proj_p) < 1e-6);

    return proj_p;
}

template class tPlane<float>;
template class tPlane<double>;

GTB_END_NAMESPACE
