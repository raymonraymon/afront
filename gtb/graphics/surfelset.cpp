
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


/*
 * Name: surfelset.cpp
 * Author: Shachar Fleishman
 * Description:
 *   
 *
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */
#include <gtb/gtb.hpp>
#ifndef WIN32
#include <gtb/error/error.hpp>
#include <gtb/graphics/surfelset.hpp>
#endif // WIN32



GTB_BEGIN_NAMESPACE

template <class T>
void tsurfel_set<T>::insert(const tsurfel_set<T>& points)
{
    //
    // Need to define a policy:
    //   1. many inserts - waste some memory
    //   2. few inserts - memory efficient.
    //
    //  The current settings is for (1).
    //
#if 0

    unsigned N = points.size();
    unsigned K = size();

    _vertices.reserve(N+K);
    if (points.has_normals()) _normals.reserve(K+N);
    if (points.has_colors()) _colors.reserve(K+N);

    for (unsigned i = 0; i < N; ++i)
    {
        _vertices.push_back(points._vertices[i]);

        if (points.has_normals()) _normals.push_back(points._normals[i]);
        if (points.has_colors()) _colors.push_back(points._colors[i]);
		if (points.has_radius()) _radius.push_back(points._radius[i]);
    }
#else
    _vertices.insert(_vertices.end(), points._vertices.begin(), points._vertices.end());
    if (points.has_normals()) _normals.insert(_normals.end(), points._normals.begin(), points._normals.end());
    if (points.has_colors()) _colors.insert(_colors.end(), points._colors.begin(), points._colors.end());
    if (points.has_radius()) _radius.insert(_radius.end(), points._radius.begin(), points._radius.end());
#endif
}

/*
 * Resize the number of vertices / normals / color to the desired size
 * 
 * The user is responsible for the results of this operation!!!
 */
template <class T>
void tsurfel_set<T>::resize(unsigned V, unsigned N, unsigned C, unsigned R)
{
    _vertices.resize(V);
    _normals.resize(N);
    _colors.resize(C);
    _radius.resize(R);
}

template <class T>
void tsurfel_set<T>::erase(unsigned first_idx, unsigned last_idx)
{
    _vertices.erase(_vertices.begin() + first_idx, _vertices.begin() + last_idx);
    if (has_normals()) _normals.erase(_normals.begin() + first_idx, _normals.begin() + last_idx);
    if (has_colors()) _colors.erase(_colors.begin() + first_idx, _colors.begin() + last_idx);
    if (has_radius()) _radius.erase(_radius.begin() + first_idx, _radius.begin() + last_idx);
}

template <class T>
void tsurfel_set<T>::erase(std::vector<unsigned>& indices)
{
    std::sort(indices.begin(), indices.end());
    int N = size();
    int K = indices.size();
    typename vertex_list::iterator lv = _vertices.end();
    typename normal_list::iterator ln; if (has_normals()) ln = _normals.end();
    typename color_list::iterator lc; if (has_colors()) lc = _colors.end();
    typename radius_list::iterator lr; if (has_radius()) lr = _radius.end();

    for (int i = K-1; i >= 0; --i)
    {
        --lv;
        std::swap(_vertices[indices[i]], *lv); 

        if (has_normals())
        {
            --ln;
            std::swap(_normals[indices[i]], *ln); 
        }
        if (has_colors())
        {
            --lc;
            std::swap(_colors[indices[i]], *lc); 
        }
        if (has_radius())
        {
            --lr;
            std::swap(_radius[indices[i]], *lr); 
        }
    }
    _vertices.erase(lv, _vertices.end());
    if (has_normals()) _normals.erase(ln, _normals.end());
    if (has_colors()) _colors.erase(lc, _colors.end());
    if (has_radius()) _radius.erase(lr, _radius.end());
}


template <class T>
void tsurfel_set<T>::assign(const tsurfel_set<T>& rhs)
{
	_vertices = rhs._vertices;
	_normals = rhs._normals;
	_colors = rhs._colors;
	_radius = rhs._radius;
}

template <class T>
void tsurfel_set<T>::compute_centroid() const
{
    tModel<T>::_centroid = tPoint3<T>::centroid(_vertices);
}

// the median of the set of points in each axis indipendently...
template <class T>
void tsurfel_set<T>::compute_median() const
{
    int N = size();
    std::vector<T> values;
    values.reserve(N);

    std::transform(_vertices.begin(), _vertices.end(), std::back_inserter(values), std::mem_fun_ref(&tPoint3<T>::x));
    std::sort(values.begin(), values.end());
    tModel<T>::_median[0] = values[N/2];

    std::transform(_vertices.begin(), _vertices.end(), values.begin(), std::mem_fun_ref(&tPoint3<T>::y));
    std::sort(values.begin(), values.end());
    tModel<T>::_median[1] = values[N/2];

    std::transform(_vertices.begin(), _vertices.end(), values.begin(), std::mem_fun_ref(&tPoint3<T>::z));
    std::sort(values.begin(), values.end());
    tModel<T>::_median[2] = values[N/2];
}

template <class T>
void tsurfel_set<T>::compute_bounding_box() const
{
    tModel<T>::_bounding_box = tBox3<T>::bounding_box(_vertices);
}

template <class T>
void tsurfel_set<T>::reserve(int count, bool rnormals, bool rcolors, bool rradius)
{
    _vertices.reserve(count);
    if (rnormals) _normals.reserve(count);
    if (rcolors) _colors.reserve(count);
    if (rradius) _radius.reserve(count);
}

template <class T>
void tsurfel_set<T>::reserve_radiuses(int count)
{
    _radius.reserve(count);
}

template <class T>
void tsurfel_set<T>::reserve_colors(int count)
{
    _colors.reserve(count);
}

template <class T>
void tsurfel_set<T>::reserve_normals(int count)
{
    _normals.reserve(count);
}

template <class T>
void tsurfel_set<T>::reserve_vertices(int count)
{
    _vertices.reserve(count);
}

template <class T>
void tsurfel_set<T>::apply(const tMatrix4<T>& M)
{
    typename vertex_list::iterator f = _vertices.begin();
    typename vertex_list::iterator l = _vertices.end();

    for (; f != l; ++f)
    {
        *f = M*(*f);
    }

    tModel<T>::invalidate_all();
}

/*
 * Convert an IndexedTriangleSet to tsurfel_set
 * if the source has no normals, new normals are NOT computed (TODO).
 */

template class tsurfel_set<float>;
template class tsurfel_set<double>;

GTB_END_NAMESPACE

