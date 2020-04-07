
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
#include <gtb/graphics/surfelset.hpp>
#include <gtb/graphics/surfelsetview.hpp>

GTB_BEGIN_NAMESPACE

unsigned surfelset_multiview::size() const
{
    unsigned K = 0;
    for (unsigned i = 0; i < _sslist.size(); ++i) K += _sslist[i]->size();
    for (unsigned i = 0; i < _ssvlist.size(); ++i) K += _ssvlist[i]->size();

    return K;
}

void surfelset_multiview::insert(surfel_set* ss)
{
    _sslist.push_back(ss);
}

void surfelset_multiview::insert(surfelset_view* ssv)
{
    _ssvlist.push_back(ssv);
}

bool surfelset_multiview::has_normals() const
{
    if (_sslist.size() > 0) return _sslist[0]->has_normals();
    else if (_ssvlist.size() > 0) return _ssvlist[0]->has_normals();
    else return false;
}

bool surfelset_multiview::has_colors() const
{
    if (_sslist.size() > 0) return _sslist[0]->has_colors();
    else if (_ssvlist.size() > 0) return _ssvlist[0]->has_colors();
    else return false;
}

bool surfelset_multiview::has_radius() const
{
    if (_sslist.size() > 0) return _sslist[0]->has_radius();
    else if (_ssvlist.size() > 0) return _ssvlist[0]->has_radius();
    else return false;
}

const Point3& surfelset_multiview::vertex(unsigned idx) const
{
    unsigned set_idx, local_idx;
    if (set_index(idx, set_idx, local_idx)) 
    {
        return _sslist[set_idx]->vertex(local_idx);
    }
    else
    {
        return _ssvlist[set_idx]->vertex(local_idx);
    }
}

const Vector3& surfelset_multiview::normal(unsigned idx) const
{
    unsigned set_idx, local_idx;
    if (set_index(idx, set_idx, local_idx)) 
    {
        return _sslist[set_idx]->normal(local_idx);
    }
    else
    {
        return _ssvlist[set_idx]->normal(local_idx);
    }
}

const ColorRgb& surfelset_multiview::vertex_color(unsigned idx) const
{
    unsigned set_idx, local_idx;
    if (set_index(idx, set_idx, local_idx)) 
    {
        return _sslist[set_idx]->vertex_color(local_idx);
    }
    else
    {
        return _ssvlist[set_idx]->vertex_color(local_idx);
    }
}

double surfelset_multiview::radius(unsigned idx) const
{
    unsigned set_idx, local_idx;
    if (set_index(idx, set_idx, local_idx)) 
    {
        return _sslist[set_idx]->radius(local_idx);
    }
    else
    {
        return _ssvlist[set_idx]->radius(local_idx);
    }
}


void surfelset_multiview::clear()
{
    _sslist.clear();
    _ssvlist.clear();
}

bool surfelset_multiview::set_index(int idx, unsigned& set_idx, unsigned& local_idx) const
{
    unsigned Kss = _sslist.size();
    set_idx = 0;
    local_idx = idx;
    while ( (set_idx < Kss) && (local_idx >= _sslist[set_idx]->size()) )
    {
        local_idx -= _sslist[set_idx]->size();
        ++set_idx;
    }
    if (set_idx < Kss) return true;

    // else
    set_idx = 0;
    unsigned Kssv = _ssvlist.size();
    while ( (set_idx < Kssv) && (local_idx >= _ssvlist[set_idx]->size()) )
    {
        local_idx -= _ssvlist[set_idx]->size();
        ++set_idx;
    }
    assert(local_idx < _ssvlist[set_idx]->size());
    return false;
}

void surfelset_multiview::compute_bounding_box() const
{
    _centroid = Point3(0,0,0);
    int N = size();
    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = vertex(i);
        for (int j = 0; j < 3; ++j)
        {
            _centroid[j] += xi[j];
        }
    }
    _centroid.scalar_scale(1.0 / N);
}

void surfelset_multiview::compute_centroid() const
{
    int N = size();
    Point3 r(-1e8,-1e8,-1e8);
    Point3 l(1e8,1e8,1e8);
    for (int i = 0; i < N; ++i)
    {
        const Point3& xi = vertex(i);
        for (int j = 0; j < 3; ++j)
        {
            l[j] = min2(l[j], xi[j]);
            r[j] = max2(r[j], xi[j]);
        }
    }
    _bounding_box = Box3(l,r);
}

GTB_END_NAMESPACE
