
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
#include <gtb/graphics/triangle3.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/triangle3.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


template <class T>
void tTriangle3<T>::subdivide(value_type max_area, vector<tTriangle3> &v) const
{
    if(area() > max_area)	{
        Point3 ab = Point3::midpoint(m_A, m_B);
        Point3 ac = Point3::midpoint(m_A, m_C);
        Point3 bc = Point3::midpoint(m_B, m_C);

        tTriangle3 t1(m_A, ab, ac);
        tTriangle3 t2(m_B, bc, ab);
        tTriangle3 t3(m_C, ac, bc);
        tTriangle3 t4(ab, bc, ac);

        t1.subdivide(max_area, v);
        t2.subdivide(max_area, v);
        t3.subdivide(max_area, v);
        t4.subdivide(max_area, v);
    } else {
        v.push_back(*this);
    }
}


template <class T>
void tTriangle3<T>::get_barycentric_coordinates(const Point3 &p,
                                                value_type &wa,
                                                value_type &wb,
                                                value_type &wc) const
{
    Vector3 ab = m_B - m_A;
    Vector3 ac = m_C - m_A;
    Vector3 ap = p - m_A;
    Vector3 bc = m_C - m_B;
    Vector3 bp = p - m_B;

    value_type den = ab.cross(ac).squared_length();
    value_type num_a = bc.cross(bp).squared_length();
    value_type num_b = ac.cross(ap).squared_length();
    value_type num_c = ab.cross(ap).squared_length();

    wa = sqrt(num_a / den);
    wb = sqrt(num_b / den);
    wc = sqrt(num_c / den);
}


// ericson closest point
template <class T>
typename tTriangle3<T>::Point3 tTriangle3<T>::closest_point(const Point3 &from) const {

	Vector3 ab = m_B - m_A;
	Vector3 ac = m_C - m_A;
	Vector3 ap = from - m_A;
	T d1 = ab.dot(ap);
	T d2 = ac.dot(ap);
	if (d1<=0 && d2<=0) return m_A;


	Vector3 bp = from - m_B;
	T d3 = ab.dot(bp);
	T d4 = ac.dot(bp);
	if (d3>=0 && d4<=d3) return m_B;


	T vc = d1*d4 - d3*d2;
	if (vc<0 && d1>=0 && d3<=0) {
		T v = d1 / (d1-d3);
		return m_A + v*ab;
	}

	Vector3 cp = from - m_C;
	T d5 = ab.dot(cp);
	T d6 = ac.dot(cp);
	if (d6>=0 && d5<=d6) return m_C;

	T vb = d5*d2 - d1*d6;
	if (vb<0 && d2>=0 && d6<=0) {
		T w = d2 / (d2-d6);
		return m_A + w*ac;
	}

	T va = d3*d6 - d5*d4;
	if (va<=0 && (d4-d3)>=0 && (d5-d6)>=0) {
		T w = (d4-d3) / ((d4-d3) + (d5-d6));
		return m_B + w*(m_C-m_B);
	}

	T denom = 1 / (va+vb+vc);
	T v = vb*denom;
	T w = vc*denom;
	return m_A + v*ab + w*ac;

}


/* Ray-Triangle Intersection Test Routines          */
/* Different optimizations of my and Ben Trumbore's */
/* code from journals of graphics tools (JGT)       */
/* http://www.acm.org/jgt/                          */
/* by Tomas Moller, May 2000                        */
template <class T>
bool tTriangle3<T>::intersect_ray(const Point3 &p0, const Vector3 &d, value_type &t) const {

	const value_type EPSILON = 0; //1e-5;
	value_type u,v;

   /* find vectors for two edges sharing vert0 */
	Vector3 edge1 = m_B - m_A;
	Vector3 edge2 = m_C - m_A;

   /* begin calculating determinant - also used to calculate U parameter */
   Vector3 pvec = d.cross(edge2);
   Vector3 qvec;

   /* if determinant is near zero, ray lies in plane of triangle */
   value_type det = edge1.dot(pvec);

   if (det > EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      Vector3 tvec  = p0 - m_A;
      
      /* calculate U parameter and test bounds */
      u = tvec.dot(pvec);
      if (u < 0 || u > det)
		 return false;
      
      /* prepare to test V parameter */
      qvec = tvec.cross(edge1);
      
      /* calculate V parameter and test bounds */
      v = d.dot(qvec);
      if (v < 0 || u+v > det)
		return false;
      
   }
   else if(det < -EPSILON)
   {
      /* calculate distance from vert0 to ray origin */
      Vector3 tvec = p0 - m_A;
      
      /* calculate U parameter and test bounds */
      u = tvec.dot(pvec);
      if (u > 0 || u < det)
		 return false;
      
      /* prepare to test V parameter */
      qvec = tvec.cross(edge1);
      
      /* calculate V parameter and test bounds */
      v = d.dot(qvec) ;
      if (v > 0 || u + v < det)
		 return false;
   }
   else return false;  /* ray is parallell to the plane of the triangle */


   value_type inv_det = 1.0 / det;

   /* calculate t, ray intersects triangle */
   t = edge2.dot(qvec) * inv_det;

   return true;
}


template <class T>
bool tTriangle3<T>::intersect_segment(const Point3 &p1, const Point3 &p2, value_type &t) const {

	Vector3 norm = normal();

	value_type dist = norm.dot(m_A);

	if ((norm.dot(p1)-dist) * (norm.dot(p2)-dist) > 0)
		return false;	// same side of plane

	Vector3 dir = p2-p1;
	value_type len = dir.length();
	dir *= 1/len;

	if (intersect_ray(p1, dir, t)) {
		t /= len;
		if (t>0 && t<1)
			return true;
	}

	return false;
}



template<class T>
bool tTriangle3<T>::Intersection(const tTriangle3<T> &tri1, const tTriangle3<T> &tri2, typename tTriangle3<T>::Point3 &i1, typename tTriangle3<T>::Point3 &i2) {

	// we want to reference the points with indices
	const typename tTriangle3<T>::Point3 tpts1[3] = { tri1.A(), tri1.B(), tri1.C() };
	const typename tTriangle3<T>::Point3 tpts2[3] = { tri2.A(), tri2.B(), tri2.C() };

	typename tTriangle3<T>::Vector3 norm1 = (tpts1[1]-tpts1[0]).cross(tpts1[2]-tpts1[0]);
	norm1.normalize();
	typename tTriangle3<T>::Vector3 norm2 = (tpts2[1]-tpts2[0]).cross(tpts2[2]-tpts2[0]);
	norm2.normalize();


	// the direction of the line of intersecting planes
	typename tTriangle3<T>::Vector3 idir = norm1.cross(norm2);
	if (idir.length()==0) return false;
	idir.normalize();

	// get a point on the line of intersection of the two planes
	typename tTriangle3<T>::Point3 ip;
	{
		typename tTriangle3<T>::Vector3 ipdir = norm1.cross(idir);

		typename tTriangle3<T>::value_type dscale = norm2.dot(ipdir);
		if (dscale==0) return false;
		ip = tpts1[0] + ((norm2.dot(tpts2[0])-norm2.dot(tpts1[0])) / dscale) * ipdir;
	}


	// see where ip + t*idir intersects with the 2 triangles
	typename tTriangle3<T>::value_type tmin=-1e34, tmax=1e34;

	typename tTriangle3<T>::Point3 tpts[2][3] = { { tpts1[0],tpts1[1],tpts1[2] }, { tpts2[0],tpts2[1],tpts2[2] } };
	typename tTriangle3<T>::Vector3 norms[2] = { norm1, norm2 };

	for (int tri=0; tri<2; tri++) {

		typename tTriangle3<T>::Vector3 udir = idir;
		typename tTriangle3<T>::Vector3 vdir = norms[tri].cross(udir);
		udir.normalize(); vdir.normalize();

		typename tTriangle3<T>::Point2 points2d[3] = { typename tTriangle3<T>::Point2(udir.dot(tpts[tri][0]-ip), vdir.dot(tpts[tri][0]-ip)),
														typename tTriangle3<T>::Point2(udir.dot(tpts[tri][1]-ip), vdir.dot(tpts[tri][1]-ip)),
														typename tTriangle3<T>::Point2(udir.dot(tpts[tri][2]-ip), vdir.dot(tpts[tri][2]-ip)) };


		bool done = false;

		for (int i=0; i<3; i++) {

			// if i is above the x axis and the other two are below..
			if ((points2d[i][1]>0 && points2d[(i+1)%3][1]<0 && points2d[(i+2)%3][1]<0) ||
				(points2d[i][1]<0 && points2d[(i+1)%3][1]>0 && points2d[(i+2)%3][1]>0)) {

				typename tTriangle3<T>::value_type alpha1 = points2d[i][1] / (points2d[i][1] - points2d[(i+1)%3][1]);
				typename tTriangle3<T>::value_type int1 = (1-alpha1)*points2d[i][0] + alpha1*points2d[(i+1)%3][0];

				typename tTriangle3<T>::value_type alpha2 = points2d[i][1] / (points2d[i][1] - points2d[(i+2)%3][1]);
				typename tTriangle3<T>::value_type int2 = (1-alpha2)*points2d[i][0] + alpha2*points2d[(i+2)%3][0];

				if (int2 < int1) {
					typename tTriangle3<T>::value_type tmp = int1;
					int1 = int2;
					int2 = tmp;
				}

				if (int1 > tmax) return false;
				if (int1 > tmin) tmin=int1;

				if (int2 < tmin) return false;
				if (int2 < tmax) tmax=int2;

				done=true;	break;
			}

		}

		if (!done)	return false;	// one of the tri's didn't intersect the other's plane

	}

	i1 = ip + tmin*idir;
	i2 = ip + tmax*idir;



	return true;
}


template class tTriangle3<float>;
template class tTriangle3<double>;

GTB_END_NAMESPACE
