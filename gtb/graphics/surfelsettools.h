
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


#ifndef __SURFELSETTOOLS_H
#define __SURFELSETTOOLS_H

#include "surfelset.hpp"
#include "surfelsetview.hpp"
#include "surfelsetmultiview.h"
#include <gtb/memory/ptrs.h>
#include <gtb/graphics/kdtree.h>

GTB_BEGIN_NAMESPACE

void convert (const IndexedTriangleSet& triangle_set, surfel_set& sset);
void point3_set2surfelset(const point3_set& p3, surfel_set& sset);

#if 0
template <class T>
tPoint3<typename T::value_type> centroid(const T& points);
#endif

//
// Helper function object to store the points of the surfelset
// in my kdtree
// Type T may be surfel_set or surfelset_view
//
template<class T>
struct GetPoint_f
{
    const T& surfels;
    GetPoint_f(const T& r_surfels) : surfels(r_surfels) {}

    const tPoint3<typename T::value_type>& operator()(unsigned idx) const { return surfels.vertex(idx); }
};

// OLD
typedef GetPoint_f<surfel_set> GetPoint_f_ss;
typedef GetPoint_f<surfelset_view> GetPoint_f_ssv;

// NEW
typedef GetPoint_f<tsurfel_set<float> > GetPoint_f_ssf;
typedef GetPoint_f<tsurfel_set<double> > GetPoint_f_ssd;
typedef GetPoint_f<tsurfelset_view<float> > GetPoint_f_ssvf;
typedef GetPoint_f<tsurfelset_view<double> > GetPoint_f_ssvd;

template <class T>
inline GetPoint_f<T> gen_GetPoint_f(const T& surfels)
{
    return GetPoint_f<T>(surfels);
}


//
// Helper function object to store point3_set in my kdtree
//
#if 0
struct GetPointp3_f
{
	const point3_set& pset;
	GetPointp3_f(const point3_set& pset_) : pset(pset_) {}

	const tPoint3& operator()(unsigned idx) const { return pset[idx]; }
};

inline GetPointp3_f gen_GetPointp3_f(const point3_set& pset)
{
    return GetPointp3_f(pset);
}
#endif

/*-------------------------------------------------------------*/

////static unsigned bpnumber = 6602;

//
// Helper function object that returns true if the 
// angle between a given surfel and the input
// surfel is not too big
// used for kdtree Traverse_if, etc.
//
struct SurfelNormalCompare
{
    const surfel_set& surfels;
    const Vector3& _n;
    double _cos_legal_angle;

    SurfelNormalCompare(const surfel_set& r_surfels, const Vector3& n, double legal_angle) : 
        surfels(r_surfels),
        _n(n),
        _cos_legal_angle(cos(legal_angle))
    {}

    bool operator()(unsigned idx) const;
};

inline SurfelNormalCompare gen_SurfelNormalCompare(
    const surfel_set& surfels, 
    const Vector3& n, 
    double legal_angle)
{
    SurfelNormalCompare snc(surfels, n, legal_angle);
    return snc;
}

//
// KDTree for surfelsets
// i.e. urfel_set or surfelset_view
//
template <class T>
class ss_kdtree
{
public:
    typedef typename T::value_type value_type;

    ss_kdtree(T& points);
    void rebuild(const tBox3<value_type>* ibbox=0);

    void extract(const tPoint3<value_type>& x, value_type radius, tsurfelset_view<value_type>& NN) const;
    void extract(const tPoint3<value_type>& x, unsigned K, tsurfelset_view<value_type>& NN) const;

    T& get_points();

    typedef KDTree<int, value_type, GetPoint_f<T> > t_surfel_tree;
    aptr<t_surfel_tree> tree;
    T& _points;
};

typedef ss_kdtree<surfel_set> kd_ss;
typedef ss_kdtree<surfelset_view> kd_ssv;

//typedef ss_kdtree<surfelset_multiview> kd_ssmv;

//
// Return the index to the point that is closest to x
//
template <class SSTYPE>
//unsigned closest_point_index(SSTYPE& cont, const Point3& x);
unsigned closest_point_index(SSTYPE& cont, const tPoint3<typename SSTYPE::value_type>& x);

GTB_END_NAMESPACE

#include "surfelsettools.inl"

#endif // __SURFELSETTOOLS_H
