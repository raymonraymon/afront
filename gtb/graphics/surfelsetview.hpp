
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


#ifndef __SURFELSETVIEW_HPP
#define __SURFELSETVIEW_HPP

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

/*
 * Name: SurfelSetView.hpp
 * Author: Shachar Fleishman
 * Description:
 *   A view of a surfelset, i.e. iterate on a subset of
 *   a surfelset.
 *
 *   The view is defined by an array of indices to the
 *   objects in the given surfel-set.
 *   I.e. the list of indices is like pointers to pointer
 *
 * TODO:
 *   create a non-const view
 *  
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */

template <class T>
class tsurfelset_view : public tModel<T>
{
public:
    typedef T value_type;

    //
    // List of indices that are visible in this view
    //
    typedef std::vector<unsigned> subset_view;
    typedef subset_view::iterator rit;
    typedef subset_view::const_iterator const_rit;

    tsurfelset_view(const tsurfelset_view& rhs);
    tsurfelset_view(const tsurfel_set<T>& surfelset);

    unsigned size() const { return _view.size(); }

    const tsurfel_set<T>& get_points() const { return _surfelset; }

    subset_view& get_view() { return _view; }
    const subset_view& get_view() const { return _view; }

    void insert(unsigned i) { 
	_view.push_back(i);  
	tModel<T>::_is_centroid_valid = false; 
    }
    void insertall();
    unsigned get_index(unsigned idx) const;
    unsigned find_index(unsigned idx) const;
    void set_index(unsigned idx, unsigned pidx);

    bool has_normals() const { return _surfelset.has_normals(); }
    bool has_colors() const { return _surfelset.has_colors(); }
    bool has_radius() const { return _surfelset.has_radius(); }

    const tPoint3<T>& vertex(rit i) const;
    const tPoint3<T>& vertex(const_rit i) const;
    const tPoint3<T>& vertex(unsigned i) const;
    const tVector3<T>& normal(unsigned i) const;
    const ColorRgb& vertex_color(unsigned i) const;
    T radius(unsigned idx) const { return _surfelset.radius(_view[idx]); }

    void clear() { 
	_view.clear(); 
	tModel<T>::_is_centroid_valid = false; 
    }
    void erase(tsurfelset_view& rhs);
    void erase(unsigned idx);
    void erase(unsigned idx1, unsigned idx2);
    void fill(int start, int end);

    typedef unsigned const_iterator;
    typedef unsigned iterator;
    iterator begin() { return 0; }
    const_iterator begin() const { return 0; }
    iterator end() { return size(); }
    const_iterator end() const { return size(); }

    rit ibegin() { return _view.begin(); }
    rit iend() { return _view.end(); }
    const_rit ibegin() const { return _view.begin(); }
    const_rit iend() const { return _view.end(); }

    tsurfelset_view& operator=(const tsurfelset_view& rhs);

    //
    // Debugging
    //
    void printme(const char* header);

protected:
    const tsurfel_set<T>& _surfelset;
    subset_view _view;

    void compute_bounding_box() const;
    void compute_centroid() const;
    void compute_median() const;
};

typedef tsurfelset_view<float> surfelset_viewf;
typedef tsurfelset_view<double> surfelset_viewd;
#if defined(REAL_IS_FLOAT)
typedef surfelset_viewf surfelset_view;
#else
typedef surfelset_viewd surfelset_view;
#endif

GTB_END_NAMESPACE

#include "surfelsetview.inl"

#endif // __SURFELSETVIEW_HPP
