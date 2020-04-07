
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


#ifndef GTB_INDEXED_TRIANGLE_SET_INCLUDED
#define GTB_INDEXED_TRIANGLE_SET_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/model.hpp>
#include <gtb/graphics/indexed_triangle.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/math/matrix4.hpp>
#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/color_rgb.hpp>
#include <gtb/graphics/triangle3.hpp>


GTB_BEGIN_NAMESPACE


class ModelRtpi;
class RtpiGrid;
class Image;
//class Camera;

template <class T>
class tIndexedTriangleSet : public tModel<T> {
public:
    typedef T value_type;
    typedef tPoint2<T> Point2;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;
    typedef tMatrix4<T> Matrix4;
    typedef tTriangle3<T> Triangle3;
    typedef tBox3<T> Box3;
    typedef tCamera<T> Camera;

    typedef std::vector<Point3> vertex_list;
    typedef std::vector<ColorRgb> color_list;
    typedef std::vector<IndexedTriangle> face_list;
    typedef std::vector<Vector3> normal_list;

    tIndexedTriangleSet();
    tIndexedTriangleSet(const tIndexedTriangleSet &its);
    explicit tIndexedTriangleSet(const ModelRtpi &rtpi_model);
    tIndexedTriangleSet(const RtpiGrid &grid, value_type max_neighbor_dist);
    virtual ~tIndexedTriangleSet();

    void clear();

    void add_vertex(const Point3 &p);
    void add_vertex(const Point3 &p, const ColorRgb &c);

    unsigned num_vertices() const;
    const Point3 &vertex(unsigned i) const;
    Point3 &vertex(unsigned i);
    const vertex_list &vertices() const;
    vertex_list &vertices();

    bool has_vertex_colors() const;
    unsigned num_vertex_colors() const;
    const ColorRgb &vertex_color(unsigned i) const;
    const color_list &vertex_colors() const;
    void need_colors();
    void need_vertex_colors() { need_colors(); }

    void need_face_colors();
    void set_face_color(unsigned idx, const ColorRgb& c);
    bool has_face_colors() const;
    const ColorRgb &face_color(unsigned i) const;
    const color_list &face_colors() const;

    unsigned num_triangles() const;
    const IndexedTriangle &indexed_triangle(unsigned i) const;
    Triangle3 triangle(unsigned i) const;
    const face_list &triangles() const;
    void flip_face(unsigned index);

    bool has_vertex_normals() const;
    unsigned num_vertex_normals() const;
    const Vector3 &vertex_normal(unsigned i) const;
    Vector3 &vertex_normal(unsigned i);
    const normal_list &vertex_normals() const;

    void read(const char *file_name);

    void read_off(FILE *fp,
                  bool subdivision_enabled = false,
                  value_type max_area = 1.0);
    void write_off_ascii(FILE *fp) const;

    void read_rtpi(FILE *fp);

    void read_binary(FILE *fp);
    void write_binary(FILE *fp) const;

    tIndexedTriangleSet *create_submodel(
        const std::vector<int> &arg_vertices,
        const std::vector<int> &triangles) const;

    void render(const Camera &cam,
                const ColorRgb &default_color /*= ColorRgb::RED */,
                value_type lod /*= 1.0*/,
                unsigned *num_points_rendered /*= NULL*/,
                bool points_as_cones /*= false*/) const;

    void render_geometry() const;

    void add(const tIndexedTriangleSet &model);

    void transform(const Matrix4 &m);

    void map_colors(const Image &img, const Camera &cam);

    void subdivide(const IndexedTriangle &t, value_type max_area);

    void insert_triangle(unsigned A, unsigned B, unsigned C,
                         bool subdivision_enabled, value_type max_area);

    void insert_triangle(const Point3 &A,
                         const Point3 &B,
                         const Point3 &C,
                         value_type max_side);

    void insert_vertex(const Point3& v);
    void set_vertex_normal(unsigned idx, const Vector3& n);
    void insert_vertex_normal(const Vector3& n);
    void insert_vertex_color(const ColorRgb& c);
    void set_vertex_color(unsigned idx, const ColorRgb& c);
    void compute_vertex_normals();
    Vector3 face_normal(unsigned face_index) const;
    bool is_valid() const;

    void scale(const Point3& origin, value_type a);
    void flip_normals();

    void set_vertex(const Point3& p, unsigned idx);

    value_type volume() const; // Compute the volume of the mesh

    // NEW: SHACHAR, cleanup
    void remove_degenerated_faces();
    void remove_double_sided_faces(bool remove_both=true);
    void remove_lonely_vertices();
    void erase_vertex(unsigned v);
    void erase_face(unsigned f);

    void edge_collaps(unsigned v1, unsigned v2);

    void insert_face(const std::vector<unsigned> &sides,
                     bool subdivison_enable,
                     value_type max_area);

    void apply(const Matrix4& M);

protected:
    void compute_bounding_box() const;
    void compute_centroid() const;

    void render_point_set(const Camera &cam,
                          const ColorRgb &default_color,
                          value_type lod,
                          unsigned *num_points_rendered,
                          bool points_as_cones) const;


    void read_off_vertex_ascii(FILE *fp, bool has_colors);
    void read_off_vertex_binary(FILE *fp, bool has_colors);
    void read_off_face_ascii(FILE *fp,
                             bool subdivision_enabled,
                             value_type max_area);
    void read_off_face_binary(FILE *fp,
                              bool subdivision_enabled,
                              value_type max_area);

    vertex_list _vertices;
    normal_list _vertex_normals;
    color_list _vertex_colors;
    color_list _face_colors;
    face_list _triangles;

public:
    static bool face_compare_less(const IndexedTriangle& ct1, const IndexedTriangle& ct2);
    static bool face_compare_eq(const IndexedTriangle& ct1, const IndexedTriangle& ct2);
};

typedef tIndexedTriangleSet<float> IndexedTriangleSetf;
typedef tIndexedTriangleSet<double> IndexedTriangleSetd;
#if defined(REAL_IS_FLOAT)
typedef IndexedTriangleSetf IndexedTriangleSet;
#else
typedef IndexedTriangleSetd IndexedTriangleSet;
#endif

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/indexed_triangle_set.ipp>
#endif

#endif // GTB_INDEXED_TRIANGLE_SET_INCLUDED
