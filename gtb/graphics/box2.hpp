
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


#ifndef GTB_BOX2_INCLUDED
#define GTB_BOX2_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/segment2.hpp>

GTB_BEGIN_NAMESPACE

template <class T>
class tBox2 {
public:
    typedef T value_type;
    enum { INSIDE, INTERSECT, OUTSIDE };
    enum Face { L, R, D, U };
    enum Vertex { LD, LU, RD, RU };

    enum VertexMasx {
	RIGHT_MASK = 1 << 1,
	UP_MASK =    1 << 0
    };


    tBox2();
    tBox2(const tBox2 &b);
    tBox2(const tPoint2<T> &min_pt, const tPoint2<T> &max_pt);
    tBox2(value_type x_min, value_type y_min,
          value_type x_max, value_type y_max);
    tBox2(const tPoint2<T> &centroid,
          value_type x_length,
          value_type y_length);
    tBox2 &operator=(const tBox2 &b);

    bool operator==(const tBox2 &b) const;
    bool operator!=(const tBox2 &b) const;
    bool is_empty() const;

    const tPoint2<T> &min_point() const;
    const tPoint2<T> &get_min_point() const;
    const tPoint2<T> &max_point() const;
    const tPoint2<T> &get_max_point() const;

    value_type x_min() const;
    value_type y_min() const;

    value_type x_max() const;
    value_type y_max() const;

    void update(const tPoint2<T>& p);
    void reset(const tPoint2<T> &min_pt, const tPoint2<T> &max_pt);
    void reset(value_type x_min, value_type y_min,
               value_type x_max, value_type y_max);
    void scale(value_type s);
    void enlarge(value_type s);

    void set_min(int axis, value_type value);
    void set_max(int axis, value_type value);

    value_type x_length() const;
    value_type y_length() const;

    value_type diagonal_length() const;
    value_type shortest_axis_length() const;
    value_type longest_axis_length() const;

    value_type area() const;

    tPoint2<T> centroid() const;
    bool contains(const tPoint2<T> &p) const;
    int classify_position(const tBox2 &b) const;
    tSegment2<T> face(unsigned i) const; // i < 4
    tLine2<T> plane(unsigned i) const; // i < 4
    tPoint2<T> vertex(unsigned i) const; // i < 4
    tPoint2<T> vertex(unsigned face, unsigned i) const; // i < 2
    tVector2<T> normal(unsigned i) const; // i < 4

    void render() const;
    void outline() const;

    void read(FILE *fp);
    void write(FILE *fp) const;

    static tBox2 bounding_box(const tPoint2<T> &a,
                              const tPoint2<T> &b,
                              const tPoint2<T> &c);
    static tBox2 bounding_box(const std::vector<tPoint2<T> > &v);
    static tBox2 make_union(const tBox2 &b1, const tBox2 &b2);

//    bool is_visible() const;

    value_type distance(const tPoint2<T>& p) const;

protected:
    bool is_order_correct() const;

    tPoint2<T> _min_pt, _max_pt;
    static const unsigned _vertex_indices[4][2];
};

GTB_GENERATE_CLASS_TYPEDEFS(Box2)

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/box2.ipp>
#endif

#endif // GTB_BOX2_INCLUDED
