
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


#ifndef _AFRONT_COMMON_H
#define _AFRONT_COMMON_H


#include <gtb/gtb.hpp>

#include <viewer/debug_draw.h>
#include <mlslib/CProjection.h>
#include <mlslib/poly2.h>


#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>

#ifdef WIN32
// for NAN
#include <xmath.h>
#define INFINITY INF
#endif

#if defined(REAL_IS_FLOAT)
	typedef float real_type;
#else
	typedef double real_type;
#endif

#ifndef WIN32
#define assertnever(s) { assert(0 && s); }
#endif

typedef gtb::tVector3<real_type> Vector3;
typedef gtb::tVector2<real_type> Vector2;
typedef gtb::tPoint3<real_type> Point3;
typedef gtb::tPoint2<real_type> Point2;
typedef gtb::tBox3<real_type> Box3;
typedef gtb::tTriangle3<real_type> Triangle3;
typedef gtb::tTriangleMesh<real_type> TriangleMesh;
typedef gtb::TriangleMeshFace TriangleMeshFace;
typedef gtb::tPlane<real_type> Plane;
typedef gtb::plane_transformation<real_type> plane_transformation;
typedef gtb::tIndexedTriangleSet<real_type> IndexedTriangleSet;
typedef gtb::AMat<real_type> AMat;
typedef gtb::AVec<real_type> AVec;
typedef gtb::tMatrix4<double> Matrix4d;

typedef gtb::tsurfel_set<real_type> surfel_set;
typedef gtb::tsurfelset_view<real_type> surfelset_view;

using gtb::afree;
using gtb::sptr;
using gtb::aptr;
using gtb::CTimerOnce;

typedef mls::Poly2<real_type> Poly2;
typedef mls::CProjection<real_type> CProjection;
typedef mls::GaussianWeightFunction<real_type> GaussianWeightFunction;

//typedef thlib::CSObject CSObject;
typedef thlib::CSObject CSObject;
typedef thlib::CSObject SLObject;


using std::vector;
using std::cerr;
//using std::cout;	// don't write to cout - it is used for camera syncronization
using std::endl;

#ifdef WIN32
typedef __int64 int64;
#define HASH_MAP64 std::map<int64,int>
#define HASH_MAP32 std::map<int,int>
#else
#include <ext/hash_map>
typedef long long int64;
class hashint64 {
    public:
    size_t operator()(int64 i) const {
	return (size_t)i;
    }
};
#define HASH_MAP64 __gnu_cxx::hash_map<int64,int, hashint64 >
#define HASH_MAP32 __gnu_cxx::hash_map<int,int>
#endif




#ifndef WIN32
#define HAS_ZLIB
#endif

#ifdef HAS_ZLIB
#include <zlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif




// sometimes we want static vectors to prevent const allocation/freeing, but want an easy way to change it since it's not thread safe
//#define static_vector(type,name) static vector<type> name; name.clear();
#define static_vector(type,name) vector<type> name; name.clear();

// wipe out an object and reconstruct it
// (call destructor and default constructor, but keep the same chunk of memory)
// good for compacting vectors
template <typename T>
void wipe_object(T &obj) {
    obj.~T();
    new (&obj) T;
}


// keep acos in it's domain
inline double acossafe(double x) {
	if (x > 1.0) return acos( 1.0);
	if (x <-1.0) return acos(-1.0);
	return acos(x);
}
inline float acossafe(float x) {
	if (x > 1.0f) return acosf( 1.0f);
	if (x <-1.0f) return acosf(-1.0f);
	return acosf(x);
}


// redraw and wait for the user to push a key
void dbgClear();
void redrawAll();
void redrawAndWait(int key=0, bool force=false);
void waitForKey(int key=0);
int GetUILastKey();

// get the last point the user selected
bool getUISelectPoint(Point3 &p, int &win);


template <typename T>
void reverse(vector<T> &v) {
	for (unsigned i=0; i<v.size()/2; i++) 
		std::swap(v[i], v[v.size()-i-1]);
}


template <typename T>
T max_element(vector<T> &v) {
    T ret = v[0];
    for (unsigned i=1; i<v.size(); i++) {
	if (v[i] > ret)
	    ret = v[i];
    }
    return ret;
}

template <typename T>
T min_element(vector<T> &v) {
    T ret = v[0];
    for (unsigned i=1; i<v.size(); i++) {
	if (v[i] < ret)
	    ret = v[i];
    }
    return ret;
}




// some basic polynomial stuff
template <typename T1, typename T2>
T2 linear(const T1 c[2], const T2 &x) {
	return ((c[0]*x) + c[1]);
}



template <typename T>
T square(const T &x) {
    return x*x;
}

template <typename T1, typename T2>
T2 quadratic(const T1 c[3], const T2 &x) {
	return (((c[0]*x) + c[1])*x + c[2]);
}



template <typename T>
T cube(const T &x) {
    return x*x*x;
}

template <typename T1, typename T2>
inline T2 cubic(const T1 c[4], const T2 &x) {
	return ((((c[0]*x) + c[1])*x + c[2])*x + c[3]);
}


template <typename T>
T clamp(T x, T low, T high) {
    if (x<low) return low;
    if (x>high) return high;
    return x;
}

template <typename T>
T pow(const T &x, int p) {
	if (p==0)
		return 1;
	
	return x*pow(x, p-1);
}



template <typename T>
bool quadratic_roots(const T c[3], T &e1, T &e2) {

	if (c[0] == 0) { // linear
		if (c[1] == 0) {
			return false; // constant
		}

		e1 = e2 = -c[2]/c[1];
		return true;
	}

	T descrim = square(c[1]) - 4*c[0]*c[2];

	if (descrim < 0)
		return false; // no real roots

	real_type t;
	if (c[1]>0) {
		t = -(c[1] + sqrt(descrim)) / 2;
	} else {
		t = -(c[1] - sqrt(descrim)) / 2;
	}

	e1 = t/c[0];
	e2 = c[2]/t;

	return true;
}


template <typename T>
bool cubic_extrema(const T c[4], T &e1, T &e2) {
	const T cderiv[3] = { 3*c[0], 2*c[1], c[2] };
	return quadratic_roots(cderiv, e1, e2);
}



// since gtb::nrran1f() is not thread safe!!
inline double myran1f(int t) {
#ifndef __linux__
  return (double)rand() / RAND_MAX;
#else

  static struct drand48_data data[100];
  static bool first = true;

  if (first) { 
    static CSObject cs;
    cs.enter();
    if (first) {
      first = false;
      for (int i=0; i<100; i++) {
		  srand48_r(i*50021, &data[i]);
      }
    }
    cs.leave();
  }


  double res;
  drand48_r(&data[t], &res);
  return res;
#endif
}


#endif // _AFRONT_COMMON_H
