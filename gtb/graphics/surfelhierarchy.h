
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


#ifndef __SURFELHIERARCHY_H
#define __SURFELHIERARCHY_H

/*
 * Name: surfelhierarchy.hpp
 * Author: Shachar Fleishman
 * Description:
 *   Hierarchy of points
 *   Basically this is a vector of surfel-sets that forms the
 *   hierarchy. The model is the union of all levels.
 *
 *   This class is mainly for debug purposes.
 *
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */

#include <vector>
#include "model.hpp"
#include <gtb/graphics/indexed_triangle_set.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/color_rgb.hpp>
#include <gtb/graphics/surfelset.hpp>
#include <gtb/memory/ptrs.h>

GTB_BEGIN_NAMESPACE

//
// An index into the hierarcy holds two components
//
struct hierarchy_index
{
public:
    hierarchy_index() {}
    hierarchy_index(unsigned a_level, unsigned a_index) :
        level(a_level), index(a_index) {}

    unsigned level;
    unsigned index;
};

class surfel_hierarchy
{
public:
    class const_iterator : public hierarchy_index
    {
    public:
        const_iterator();
        const_iterator(const const_iterator& rhs);
        const_iterator(const surfel_hierarchy* hierarchy, unsigned a_level, unsigned a_index);

        const_iterator& operator++();
        bool operator==(const const_iterator& rhs);
        bool operator!=(const const_iterator& rhs);

    private:
        const surfel_hierarchy* _p_hierarchy;
    };

    const Point3& vertex(const hierarchy_index& hi) const;
    const Point3& vertex(unsigned level, unsigned index) const;
    const Point3& vertex(const_iterator& i) const;
    hierarchy_index insert_vertex(unsigned level, const Point3& v);

    const Vector3& normal(const hierarchy_index& hi) const;
    const Vector3& normal(unsigned level, unsigned index) const;
    const Vector3& normal(const_iterator& i) const;
    void insert_normal(unsigned level, const Vector3& n);

    const ColorRgb& vertex_color(const hierarchy_index& hi) const;

    void new_level();
    unsigned num_levels() const { return _hierarchy.size(); }
    const surfel_set& level(unsigned i) const { return *_hierarchy[i]; }

	unsigned size_() const; // for debug only
	void clear_level(unsigned level_);

    const_iterator begin() const;
    const_iterator end() const;

    void render() const;
    void render(unsigned level) const;
    void render(unsigned level1, unsigned level2) const;

    bool has_normals() const { return _hierarchy[0]->has_normals(); }
    bool has_colors() const { return _hierarchy[0]->has_colors(); }

	//
	// IO
	//
	void write_as_obj(const char* name) const;

private:
    typedef std::vector<auto_ptr2<surfel_set> > t_hierarchy;
    t_hierarchy _hierarchy;
};

class surfel_hierarchy_view
{
public:
    //
    // List of indices that are visible in this view
    //
    typedef std::vector<hierarchy_index> subset_view;

    explicit surfel_hierarchy_view(const surfel_hierarchy& hierarchy);

    void insert(const hierarchy_index& hi) { _view.push_back(hi); }
    void insert(unsigned level, unsigned index) { _view.push_back(hierarchy_index(level, index)); }
    const hierarchy_index& get_index(unsigned i) const { return _view[i]; }

    const Point3& vertex(unsigned i) const;
    const Vector3& normal(unsigned i) const;
    const ColorRgb& vertex_color(unsigned i) const;

    unsigned size() const { return _view.size(); }

    void clear() { _view.clear(); }

    subset_view& get_view() { return _view; }
    const surfel_hierarchy& get_hierarchy() { return _hierarchy; }

    bool has_normals() const { return _hierarchy.has_normals(); }
    bool has_colors() const { return _hierarchy.has_colors(); }

protected:
    const surfel_hierarchy& _hierarchy;
    subset_view _view;
};


//
// Helper function object to store the points of the surfelset
// in my kdtree
//
struct HierarchyGetPoint_f
{
    const surfel_hierarchy& _hierarchy;
    HierarchyGetPoint_f(const surfel_hierarchy& hierarchy) :
        _hierarchy(hierarchy) {}

    const Point3& operator()(const hierarchy_index& idx) const;
};

inline HierarchyGetPoint_f gen_HierarchyGetPoint_f(const surfel_hierarchy& hierarchy)
{
    return HierarchyGetPoint_f(hierarchy);
}

//
// Helper function object that returns true if the 
// angle between a given surfel and the input
// surfel is not too big
// used for kdtree Traverse_if, etc.
//
struct HierarchyNormalCompare
{
    const surfel_hierarchy& _hierarchy;
    const Vector3& _n;
    double _cos_legal_angle;

    HierarchyNormalCompare(const surfel_hierarchy& hierarchy, const Vector3& n, double legal_angle) : 
        _hierarchy(hierarchy),
        _n(n),
        _cos_legal_angle(cos(legal_angle))
    {}

    bool operator()(const hierarchy_index& idx) const;
};

inline HierarchyNormalCompare gen_HierarchyNormalCompare(
    const surfel_hierarchy& hierarchy,
    const Vector3& n, 
    double legal_angle)
{
    HierarchyNormalCompare hnc(hierarchy, n, legal_angle);
    return hnc;
}

#include "surfelhierarchy.inl"

GTB_END_NAMESPACE

#endif // __SURFELHIERARCHY_H
