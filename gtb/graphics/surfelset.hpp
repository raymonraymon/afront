
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


#ifndef __SURFELSET_HPP
#define __SURFELSET_HPP

/*
 * Name: SurfelSet.hpp
 * Author: Shachar Fleishman
 * Description:
 *   Surfel set class: a set of points with attributes
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

GTB_BEGIN_NAMESPACE

template<class T>
class tsurfel_set : public tModel<T>
{
public:
    typedef T value_type;

    typedef tBox3<T> Box3;

    typedef std::vector<tPoint3<T> > vertex_list;
    typedef std::vector<tVector3<T> > normal_list;
    typedef std::vector<ColorRgb> color_list;
    typedef std::vector<T> radius_list;

    tsurfel_set() {}

    virtual ~tsurfel_set() {}

    void resize(unsigned V, unsigned N=0, unsigned C=0, unsigned R=0);

    unsigned size() const;
    void clear();

    tPoint3<T>& vertex(unsigned idx);
    const tPoint3<T>& vertex(unsigned idx) const;
    const vertex_list &vertices() const;
    vertex_list &vertices();
    void set_vertex(unsigned i, const tPoint3<T>& v);
    void insert_vertex(const tPoint3<T>& v);
    void insert_vertex(const tPoint3<T>& v, const tVector3<T>& n);
    void reserve_vertices(int count);

    bool has_normals() const;
    const tVector3<T> &normal(unsigned i) const;
    tVector3<T> &normal(unsigned i);
    const normal_list &normals() const;
    normal_list &normals();
    void set_normal(unsigned i, const tVector3<T>& n);
    void insert_normal(const tVector3<T>& n);
    void clear_normals();
    void reserve_normals(int count);

    bool has_colors() const;
    const ColorRgb &vertex_color(unsigned i) const;
    ColorRgb &vertex_color(unsigned i);
    const color_list& vertex_colors() const;
    color_list& vertex_colors();
    void set_color(unsigned i, const ColorRgb& c);
    void insert_color(const ColorRgb& c);
    void allocate_colors();
    void clear_colors();
    void reserve_colors(int count);

    bool has_radius() const;
    void insert_radius(T r);
    void set_radius(unsigned idx, T r);
    T radius(unsigned idx) const;
    void clear_radius();
    void reserve_radiuses(int count);
    radius_list& radiuses();
    const radius_list& radiuses() const;

    void reserve(int count, bool rnormals=false, bool rcolors=false, bool rradius=false);

    void insert(const tsurfel_set& points);

    void erase(unsigned idx);
    void erase(unsigned first_idx, unsigned last_idx);
    void erase(std::vector<unsigned>& indices);

    void assign(const tsurfel_set& rhs);

    void apply(const tMatrix4<T>& M);

    void compute_bounding_box() const;
    void compute_centroid() const;
    void compute_median() const;

protected:
    void update_bounding_box(const tPoint3<T>& v);
    void update_centroid(const tPoint3<T>& v);
    
    vertex_list _vertices;
    normal_list _normals;
    color_list _colors;
    radius_list _radius;
    
};

typedef tsurfel_set<float> surfel_setf;
typedef tsurfel_set<double> surfel_setd;
#if defined(REAL_IS_FLOAT)
typedef surfel_setf surfel_set;
#else
typedef surfel_setd surfel_set;
#endif

typedef std::vector<Point2> point2_set;
typedef std::vector<Point3> point3_set;

GTB_END_NAMESPACE

#include "surfelset.inl"

#endif // __SURFELSET_HPP
