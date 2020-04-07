
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
#include "surfelsettools.h"

GTB_BEGIN_NAMESPACE

void convert (const IndexedTriangleSet& triangle_set, surfel_set& sset)
{
	typedef IndexedTriangleSet::vertex_list::const_iterator cvi;
	typedef IndexedTriangleSet::normal_list::const_iterator cni;

	int n_vertices = triangle_set.num_vertices();
	int n_normals = triangle_set.num_vertex_normals();
	bool has_normals = n_vertices == n_normals;

	if (!has_normals)
#ifndef NO_EXCEPTIONS
	  throw CErr("convert, IndexedTriangleSet with not vertex normals");
#else
	  return;
#endif

	cvi fv = triangle_set.vertices().begin();
	cvi lv = triangle_set.vertices().end();
	cni fn = triangle_set.vertex_normals().begin();

	for (; fv != lv; ++fv, ++fn)
	{
            sset.insert_vertex(*fv, *fn);
	}
}


/*
 * inserts a point3_set to a surfelset
 */
void point3_set2surfelset(const point3_set& p3, surfel_set& sset)
{
    point3_set::const_iterator f = p3.begin();
    point3_set::const_iterator l = p3.end();
    for (; f != l; ++f)
    {
        sset.insert_vertex(*f);
    }
}



bool SurfelNormalCompare::operator()(unsigned idx) const
{
////    if (idx == bpnumber) DebugBreak(); // DEBUG ONLY
    double cos_angle = _n.dot(surfels.normal(idx));
    if (cos_angle > _cos_legal_angle) return true;
    else return false;
};



//
// Instantiate for the required types
//
/*
template class ss_kdtree<surfel_setf>;
template class ss_kdtree<surfel_setd>;
template class ss_kdtree<surfelset_viewf>;
template class ss_kdtree<surfelset_viewd>;
*/
//template class ss_kdtree<surfelset_multiview>;

/*----------------------------------------------*/
template <class T>
tPoint3<typename T::value_type> centroid(const T& points)
//Point3 centroid(const T& points)
{
    int N = points.size();
    tPoint3<typename T::value_type> center(0,0,0);
    for (int i = 0; i < N; ++i)
    {
        center.add(points.vertex(i));
    }
    center.scalar_scale(1.0/N);
    return center;
}

#if 0
template tPoint3<float> centroid(const surfel_setf& points);
template tPoint3<double> centroid(const surfel_setd& points);
template tPoint3<float> centroid(const surfelset_viewf& points);
template tPoint3<double> centroid(const surfelset_viewd& points);
#else
// Instantiation hack to overcome a bug in vc build 13.10.3077
void centroid_instantiator()
{
    surfel_setf sf;
    surfel_setd sd;
    surfelset_viewf ssf(sf);
    surfelset_viewd ssd(sd);
    centroid(sf);
    centroid(sd);
    centroid(ssf);
    centroid(ssd);
}
#endif

/*----------------------------------------------*/
//
// Return the index to the point that is closest to x
//
template <class SSTYPE>
unsigned closest_point_index(SSTYPE& cont, const tPoint3<typename SSTYPE::value_type>& x)
{
    unsigned N = cont.size();
    if (N == 0) return -1;

    unsigned closest_idx = 0;
    typename SSTYPE::value_type mind2 = (cont.vertex(0u)-x).squared_length();
    for (unsigned i = 1; i < N; ++i)
    {
        double d2 = (cont.vertex(i) - x).squared_length();
        if (d2 < mind2)
        {
            mind2 = d2;
            closest_idx = i;
        }
    }
    return closest_idx;
}

//template unsigned closest_point_index<surfel_set>(surfel_set& cont, const Point3& x);
//template unsigned closest_point_index<surfelset_view>(surfelset_view& cont, const Point3& x);
#if 0
template unsigned closest_point_index<surfelset_viewf>(surfel_setf& cont, const Point3f& x);
template unsigned closest_point_index<surfelset_viewd>(surfel_setd& cont, const Point3d& x);
template unsigned closest_point_index<surfelset_viewf>(surfelset_viewf& cont, const Point3f& x);
template unsigned closest_point_index<surfelset_viewd>(surfelset_viewd& cont, const Point3d& x);
#else
// Instantiation hack to overcome a bug in vc build 13.10.3077
void cpi_instantiator()
{
    surfel_setf sf;
    surfel_setd sd;
    surfelset_viewf ssf(sf);
    surfelset_viewd ssd(sd);
    Point3f pf;
    Point3d pd;
    unsigned x;
    x = closest_point_index(sf, pf);
    x = closest_point_index(sd, pd);
    x = closest_point_index(ssf, pf);
    x = closest_point_index(ssd, pd);
}

#endif

GTB_END_NAMESPACE
