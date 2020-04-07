
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


#ifndef __SURFELSETMULTIVIEW_H
#define __SURFELSETMULTIVIEW_H

#include <gtb/common.hpp>
#include "surfelset.hpp"
#include "surfelsetview.hpp"

GTB_BEGIN_NAMESPACE

/*
 * Name: SurfelSetMultiView.h
 * Author: Shachar Fleishman
 * Description:
 *   Like SurfelSetView, but for multiple sets
 *
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */

class surfelset_multiview : public Model
{
public:
    surfelset_multiview() : Model() {}

    //
    // List of indices that are visible in this view
    //
    typedef std::vector<surfel_set*> t_sslist;
    typedef std::vector<surfelset_view*> t_ssvlist;

    unsigned size() const;

    void insert(surfel_set* ss);
    void insert(surfelset_view* ssv);

    bool has_normals() const;
    bool has_colors() const;
    bool has_radius() const;

    const Point3& vertex(unsigned idx) const;
    const Vector3& normal(unsigned idx) const;
    const ColorRgb& vertex_color(unsigned idx) const;
    double radius(unsigned idx) const;

    void clear();

protected:
    t_sslist _sslist;
    t_ssvlist _ssvlist;

    void compute_bounding_box() const;
    void compute_centroid() const;

    /*
     * Give an index, compute the index to the correct subset
     * input:
     *   idx
     *
     * Output
     *    set_idx - the index in _ssvlist or _sslist
     *    local_idx - the index in _sslist[set_idx] or _ssvlist
     * Return:
     *    true - the set_idx is to _sslist
     *    false - the set_idx is to _ssvlist
     */
    bool set_index(int idx, unsigned& set_idx, unsigned& local_idx) const;
};

GTB_END_NAMESPACE

#endif // __SURFELSETMULTIVIEW_H
