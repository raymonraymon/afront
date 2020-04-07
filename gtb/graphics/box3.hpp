
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


#ifndef GTB_BOX3_INCLUDED
#define GTB_BOX3_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/polygon3.hpp>

GTB_BEGIN_NAMESPACE


//class Vector3;
//class Plane;
//class Polygon3;

template <class T> class tPolygon3;
template <class T> class tPlane;

template<class T>
class tBox3 {
public:
    typedef T value_type;
    enum { INSIDE, INTERSECT, OUTSIDE };

    enum Face { L, R, D, U, B, F };

    enum Vertex { LDB,
                  LDF,
                  LUB,
                  LUF,
                  RDB,
                  RDF,
                  RUB,
                  RUF
    };

    enum VertexMask {
        RIGHT_MASK =	1 << 2,
        UP_MASK =	1 << 1,
        FRONT_MASK =	1 << 0
    };

    tBox3();
    tBox3(const tBox3 &b);
    tBox3(const tPoint3<T> &min_pt, const tPoint3<T> &max_pt);
    tBox3(value_type x_min, value_type y_min, value_type z_min,
          value_type x_max, value_type y_max, value_type z_max);
    tBox3(const tPoint3<T> &centroid,
          value_type x_length,
          value_type y_length, 
          value_type z_length);
    tBox3 &operator=(const tBox3 &b);

    bool operator==(const tBox3 &b) const;
    bool operator!=(const tBox3 &b) const;
    bool is_empty() const;

    const tPoint3<T> &min_point() const;
    const tPoint3<T> &get_min_point() const;
    const tPoint3<T> &max_point() const;
    const tPoint3<T> &get_max_point() const;

    value_type x_min() const;
    value_type y_min() const;
    value_type z_min() const;

    value_type x_max() const;
    value_type y_max() const;
    value_type z_max() const;

    void update(const tPoint3<T>& p);
    void reset(const tPoint3<T> &min_pt, const tPoint3<T> &max_pt);
    void reset(value_type x_min, value_type y_min, value_type z_min,
               value_type x_max, value_type y_max, value_type z_max);
    void scale(value_type s);
    void enlarge(value_type s);

    void set_min(int axis, value_type value);
    void set_max(int axis, value_type value);

    value_type x_length() const;
    value_type y_length() const;
    value_type z_length() const;

    value_type diagonal_length() const;
    value_type shortest_axis_length() const;
    value_type longest_axis_length() const;

    value_type volume() const;

    tPoint3<T> centroid() const;
    bool contains(const tPoint3<T> &p) const;
    int classify_position(const tBox3 &b) const;
    tPolygon3<T> face(unsigned i) const; // i < 6
    tPlane<T> plane(unsigned i) const; // i < 6
    tPoint3<T> vertex(unsigned i) const; // i < 8
    tPoint3<T> vertex(unsigned face, unsigned i) const; // i < 4
    tVector3<T> normal(unsigned i) const; // i < 6

    void render() const;
    void outline() const;

    void read(FILE *fp);
    void write(FILE *fp) const;

    static tBox3 bounding_box(const tPoint3<T> &a,
                              const tPoint3<T> &b,
                              const tPoint3<T> &c);
    static tBox3 bounding_box(const std::vector<tPoint3<T> > &v);
    static tBox3 make_union(const tBox3 &b1, const tBox3 &b2);
    static tBox3 make_intersection(const tBox3 &b1, const tBox3 &t2);

//    bool is_visible() const;

    value_type distance(const tPoint3<T>& p) const;

protected:
    bool is_order_correct() const;

    tPoint3<T> _min_pt, _max_pt;
    static const unsigned _vertex_indices[6][4];
};

GTB_GENERATE_CLASS_TYPEDEFS(Box3)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/box3.ipp>
#endif

#endif // GTB_BOX3_INCLUDED
