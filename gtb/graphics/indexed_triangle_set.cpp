
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
#ifndef WIN32
#include <gtb/graphics/indexed_triangle_set.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/cone.hpp>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/viewport.hpp>
#include <gtb/graphics/view.hpp>
#include <gtb/graphics/image.hpp>
#include <gtb/graphics/material.hpp>
#include <gtb/graphics/model_rtpi.hpp>
#include <gtb/graphics/rtpi_grid.hpp>
#include <gtb/graphics/triangle3.hpp>
#include <gtb/io/io.hpp>
#include <gtb/error/error.hpp>
#include <map>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/indexed_triangle_set.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE

using namespace std;

template <class T>
tIndexedTriangleSet<T>::tIndexedTriangleSet()
{
}


template <class T>
tIndexedTriangleSet<T>::tIndexedTriangleSet(const tIndexedTriangleSet &its)
    : tModel<T>(),
      _vertices(its._vertices),
      _vertex_normals(its._vertex_normals),
      _vertex_colors(its._vertex_colors),
      _face_colors(its._face_colors),
      _triangles(its._triangles)
{
}


template <class T>
tIndexedTriangleSet<T>::tIndexedTriangleSet(const ModelRtpi &model)
{
    int intensity_range = model.max_intensity() - model.min_intensity();

    // for each strip
    for (unsigned i = 0; i < model.num_strips(); i++) {
        const RtpiStrip &strip = model.strip(i);

        // for each point
        for (unsigned j = 0; j < strip.num_points(); j++) {
            const Rtpi &p = strip[j];

            // insert point
            _vertices.push_back(Point3(p.point()));

            // insert color
            value_type c = ((value_type) (p.i() - model.min_intensity()) /
                            (value_type) intensity_range);
            _vertex_colors.push_back(ColorRgb(c, c, c));
        }
    }
}


template <class T>
void tIndexedTriangleSet<T>::read_rtpi(FILE *fp)
{
    assert(fp != NULL);

    // read header
    unsigned n_strips, dummy;
    ModelRtpi::read_header(fp, n_strips, dummy);

    vector<int> intensities;

    // for each strip
    for (unsigned i = 0; i < n_strips; i++) {
        unsigned n;
        read_unsigned(&n, fp);

        // for each point
        for (unsigned j = 0; j < n; j++) {
            Rtpi p;
            p.read(fp);

            // insert point
            _vertices.push_back(Point3(p.point()));

            // save intensity
            intensities.push_back(p.i());
        }
    }

    // find intensity range
    int max_intensity = 0;
    int min_intensity = INT_MAX;
    for (unsigned i = 0; i < intensities.size(); i++) {
        int intensity = intensities[i];
        if (intensity < min_intensity) {
            min_intensity = intensity;
        } else if (intensity > max_intensity) {
            max_intensity = intensity;
        }
    }

    // assign colors
    int intensity_range = max_intensity - min_intensity;
    assert(_vertices.size() == intensities.size());
    for (unsigned i = 0; i < _vertices.size(); i++) {
        int intensity = intensities[i];
        value_type c = ((value_type) (intensity - min_intensity) /
                        (value_type) intensity_range);
        _vertex_colors.push_back(ColorRgb(c, c, c));
    }
}


template <class T>
void tIndexedTriangleSet<T>::insert_triangle(const Point3 &A,
                                             const Point3 &B,
                                             const Point3 &C,
                                             value_type max_side)
{
    if (Point3::distance(A, B) <= max_side &&
        Point3::distance(B, C) <= max_side &&
        Point3::distance(A, C) <= max_side) {
        unsigned n = _vertices.size();
        IndexedTriangle t(n, n + 1, n + 2);
        _triangles.push_back(t);
        _vertices.push_back(A);
        _vertices.push_back(B);
        _vertices.push_back(C);
    }
}


template <class T>
void tIndexedTriangleSet<T>::insert_triangle(
    unsigned A,
    unsigned B,
    unsigned C,
    bool subdivision_enabled,
    value_type max_area)
{
    IndexedTriangle it(A, B, C);
    if (subdivision_enabled) {
        subdivide(it, max_area);
    } else {
        _triangles.push_back(it);
    }
}


template <class T>
void tIndexedTriangleSet<T>::insert_face(
    const vector<unsigned> &sides,
    bool subdivision_enabled,
    value_type max_area)
{
    unsigned n_sides = sides.size();
    bool has_colors = _vertex_colors.size() > 0;

    if (n_sides == 3) {
        insert_triangle(sides[0], sides[1], sides[2],
                        subdivision_enabled, max_area);
    } else {
#if 0
        // find centroid
        vector<Point3> polygon(n_sides);
        for (unsigned j = 0; j < n_sides; j++) {
            polygon[j] = _vertices[sides[j]];
        }
        Point3 pc = Point3::centroid(polygon);
        int ipc = _vertices.size();
        _vertices.push_back(pc);

        // compute color for centroid
        if (has_colors)
        {
            ColorRgb crgb;
            for (unsigned j = 0; j < n_sides; j++) {
                crgb += _vertex_colors[sides[j]] / n_sides;
            }
            _vertex_colors.push_back(crgb);
        }

        // split face into triangles
        for (unsigned j = 0; j < n_sides; j++) {
            unsigned k = (j + 1) % n_sides;
            insert_triangle(ipc, sides[j], sides[k],
                            subdivision_enabled, max_area);
        }
#endif
        int new_vertices_indices[3];
        for (int i = 0; i < 3; ++i)
        {
            int i1 = (i+1)%3;
            Point3 p = Point3::midpoint(_vertices[sides[i]], _vertices[sides[i1]]);
            new_vertices_indices[i] = _vertices.size();
            _vertices.push_back(p);
        }
        insert_triangle(sides[0], new_vertices_indices[0], new_vertices_indices[2], subdivision_enabled, max_area);
        insert_triangle(sides[1], new_vertices_indices[1], new_vertices_indices[0], subdivision_enabled, max_area);
        insert_triangle(sides[2], new_vertices_indices[2], new_vertices_indices[1], subdivision_enabled, max_area);
        insert_triangle(new_vertices_indices[0], new_vertices_indices[1], new_vertices_indices[2], subdivision_enabled, max_area);
    }
}


template <class T>
tIndexedTriangleSet<T>::tIndexedTriangleSet(const RtpiGrid &grid,
                                            value_type max_neighbor_dist)
{
    for (unsigned k = 0; k < grid.num_rows() - 1; k++) {
        for (unsigned l = 0; l < grid.num_cols() - 1; l++) {
            insert_triangle(Point3(grid.sample(k, l).point()),
                            Point3(grid.sample(k + 1, l).point()),
                            Point3(grid.sample(k + 1, l + 1).point()),
                            max_neighbor_dist);

            insert_triangle(Point3(grid.sample(k, l).point()),
                            Point3(grid.sample(k, l + 1).point()),
                            Point3(grid.sample(k + 1, l + 1).point()),
                            max_neighbor_dist);
        }
    }
}


template <class T>
tIndexedTriangleSet<T>::~tIndexedTriangleSet()
{
}


template <class T>
void tIndexedTriangleSet<T>::read_binary(FILE *fp)
{
    assert(fp != NULL);

    tModel<T>::_bounding_box.read(fp);
    tModel<T>::_centroid.read(fp);

    // read header
    unsigned n_vertices;
    bool has_colors;
    unsigned n_triangles;
    bool has_normals;

    read_unsigned(&n_vertices, fp);
    read_bool(&has_colors, fp);
    read_unsigned(&n_triangles, fp);
    read_bool(&has_normals, fp);

    // read vertices, colors, and normals
    _vertices.resize(n_vertices);
    if (has_colors) {
        _vertex_colors.resize(n_vertices);
    }
    if (has_normals) {
        _vertex_normals.resize(n_vertices);
    }
    for (unsigned i = 0; i < n_vertices; i++) {
        _vertices[i].read(fp);
        if (has_colors) {
            _vertex_colors[i].read(fp);
        }
        if (has_normals) {
            _vertex_normals[i].read(fp);
        }
    }

    // read triangles
    _triangles.resize(n_triangles);
    for (unsigned i = 0; i < n_triangles; i++) {
        _triangles[i].read(fp);
    }
}


template <class T>
void tIndexedTriangleSet<T>::write_binary(FILE *fp) const
{
    assert(fp != NULL);

    tModel<T>::_bounding_box.write(fp);
    tModel<T>::_centroid.write(fp);

    // write header
    unsigned n_vertices = _vertices.size();
    bool has_colors = _vertex_colors.size() > 0;
    unsigned n_triangles = _triangles.size();
    bool has_normals = _vertex_normals.size() > 0;

    write_unsigned(n_vertices, fp);
    write_bool(has_colors, fp);
    write_unsigned(n_triangles, fp);
    write_bool(has_normals, fp);

    assert(has_colors ? n_vertices == _vertex_colors.size() : 1);
    assert(has_normals ? n_vertices == _vertex_normals.size() : 1);

    // write vertices, colors, and normals
    for (unsigned i = 0; i < n_vertices; i++) {
        _vertices[i].write(fp);
        if (has_colors) {
            _vertex_colors[i].write(fp);
        }
        if (has_normals) {
            _vertex_normals[i].write(fp);
        }
    }

    // write triangles
    for (unsigned i = 0; i < n_triangles; i++) {
        _triangles[i].write(fp);
    }
}


template <class T>
tIndexedTriangleSet<T> *tIndexedTriangleSet<T>::create_submodel(
    const vector<int> &arg_vertices,
    const vector<int> &arg_triangles) const
{
    tIndexedTriangleSet *model = new tIndexedTriangleSet();
    assert(model != NULL);

    bool has_colors = _vertex_colors.size() > 0;
    bool has_normals = _vertex_normals.size() > 0;

    if (arg_triangles.size() == 0) {
        for (unsigned i = 0; i < arg_vertices.size(); i++) {
            int j = arg_vertices[i];
            assert(j >= 0);
            assert(j < (int) _vertices.size());
            model->_vertices.push_back(_vertices[j]);
            if (has_colors) {
                model->_vertex_colors.push_back(
                    _vertex_colors[j]);
            }
            if (has_normals) {
                model->_vertex_normals.push_back(
                    _vertex_normals[j]);
            }
        }
    } else {
        for (unsigned i = 0; i < arg_triangles.size(); i++) {
            int j = arg_triangles[i];
            assert(j >= 0);
            assert(j < (int) _triangles.size());
            unsigned n = model->_vertices.size();
            model->_triangles.push_back(
                IndexedTriangle(n, n + 1, n + 2));

            const IndexedTriangle &t = _triangles[j];
            model->_vertices.push_back(_vertices[t.A()]);
            model->_vertices.push_back(_vertices[t.B()]);
            model->_vertices.push_back(_vertices[t.C()]);

            if (has_colors) {
                model->_vertex_colors.push_back(
                    _vertex_colors[t.A()]);
                model->_vertex_colors.push_back(
                    _vertex_colors[t.B()]);
                model->_vertex_colors.push_back(
                    _vertex_colors[t.C()]);
            }

            if (has_normals) {
                model->_vertex_normals.push_back(
                    _vertex_normals[t.A()]);
                model->_vertex_normals.push_back(
                    _vertex_normals[t.B()]);
                model->_vertex_normals.push_back(
                    _vertex_normals[t.C()]);
            }
        }
    }
    return model;
}


template <class T>
void tIndexedTriangleSet<T>::render_point_set(const Camera &cam,
                                              const ColorRgb &default_color,
                                              value_type lod,
                                              unsigned *num_points_rendered,
                                              bool points_as_cones) const
{
    if (0.0 == lod) {
        return;
    }
    if (num_vertex_colors() == 0) {
        default_color.load();
    }
    if (points_as_cones) {
        // FIXME
        value_type cone_height = 1.2;
        value_type cone_radius = 0.2;
        unsigned j = (unsigned) (1.0 / lod);
        Vector3 towards = cam.towards();
        for (unsigned i = 0; i < num_vertices(); i += j) {
            if (NULL != num_points_rendered) {
                (*num_points_rendered)++;
            }
            if (num_vertex_colors() > 0) {
                vertex_color(i).load();
            }
#if 0
            Cone c(vertex(i),
                   vertex(i) + cone_height * towards,
                   cone_radius);
            c.render();
#endif
        }
        return;
    }
    if (has_vertex_colors() && (!has_vertex_normals())) {
        unsigned j = (unsigned) (1.0 / lod);
        unsigned n = 0;
        glBegin(GL_POINTS);
        for (unsigned i = 0; i < num_vertices(); i += j) {
            n++;
            vertex_color(i).load();
            vertex(i).load();
        }
        glEnd();
        if (NULL != num_points_rendered) {
            *num_points_rendered = n;
        }
    } else {
        glBegin(GL_POINTS);
        unsigned j = (unsigned) (1.0 / lod);
        for (unsigned i = 0; i < num_vertices(); i += j) {
            if (NULL != num_points_rendered) {
                (*num_points_rendered)++;
            }
            if (has_vertex_colors()) {
                vertex_color(i).load();
            }
            if (has_vertex_normals()) {
                vertex_normal(i).load_as_normal();
            }
            vertex(i).load();
        }
        glEnd();
    }
}


template <class T>
void tIndexedTriangleSet<T>::render(const Camera &cam,
                                    const ColorRgb &default_color,
                                    value_type lod,
                                    unsigned *num_points_rendered,
                                    bool points_as_cones) const
{
    if (num_triangles() == 0) {
        render_point_set(cam, default_color, lod, num_points_rendered,
                         points_as_cones);
        return;
    }

    if (!has_vertex_normals()) {
        GTB_ERROR("not implemented yet");
    }

    if (num_vertex_colors() == 0) {
        default_color.load();
        glBegin(GL_TRIANGLES);
        for (unsigned i = 0; i < num_triangles(); i++) {
            const IndexedTriangle &ti = indexed_triangle(i);

            vertex_normal(ti[0]).load_as_normal();
            vertex(ti[0]).load();

            vertex_normal(ti[1]).load_as_normal();
            vertex(ti[1]).load();

            vertex_normal(ti[2]).load_as_normal();
            vertex(ti[2]).load();
        }
        glEnd();
    } else {
        glBegin(GL_TRIANGLES);
        for (unsigned i = 0; i < num_triangles(); i++) {
            const IndexedTriangle &ti = indexed_triangle(i);

            vertex_color(ti[0]).load();
            vertex_normal(ti[0]).load_as_normal();
            vertex(ti[0]).load();

            vertex_color(ti[1]).load();
            vertex_normal(ti[1]).load_as_normal();
            vertex(ti[1]).load();

            vertex_color(ti[2]).load();
            vertex_normal(ti[2]).load_as_normal();
            vertex(ti[2]).load();
        }
        glEnd();
    }
}


template <class T>
void tIndexedTriangleSet<T>::render_geometry() const
{
    if (num_triangles() == 0) {
        glBegin(GL_POINTS);
        for (unsigned i = 0; i < num_vertices(); i ++) {
            vertex(i).load();
        }
        glEnd();
    } else {
        glBegin(GL_TRIANGLES);
        for (unsigned i = 0; i < num_triangles(); i++) {
            const IndexedTriangle &ti = indexed_triangle(i);
            Vector3 n = (vertex(ti[1])- vertex(ti[0])).cross(vertex(ti[2]) - vertex(ti[0]));
            n.load_as_normal();
            vertex(ti[0]).load();
            vertex(ti[1]).load();
            vertex(ti[2]).load();
        }
        glEnd();
    }
}


template <class T>
void tIndexedTriangleSet<T>::compute_bounding_box() const
{
    tModel<T>::_bounding_box = Box3::bounding_box(_vertices);
}


template <class T>
void tIndexedTriangleSet<T>::compute_centroid() const
{
    tModel<T>::_centroid = Point3::centroid(_vertices);
}


static bool is_empty_line(const char* s)
{
    assert(s != 0);
    const char* p = s;
    while ((*p != '\0') && isspace(*p)) ++p;
    if (*p == '\0') return true;
    else return false;
}

static void read_off_header(
    FILE *fp,
    bool &has_colors,
    unsigned &n_vertices,
    unsigned &n_faces,
    unsigned &n_edges,
    bool &binary)
{
    assert(fp);
    int c;
    while (((c = getc(fp)) != EOF) && (c != 'C') && (c != 'O')) {
        ;
    }
    assert((c == 'C') || (c == 'O'));
    if (c == 'C') {
        has_colors = true;
        c = getc(fp);
    }
    assert(c == 'O');
    c = getc(fp);
    assert(c == 'F');
    c = getc(fp);
    assert(c == 'F');
    while (((c = getc(fp)) != EOF) && (c != '\n')) {
        if (c == 'B') {
            c = getc(fp);
            assert(c == 'I');
            c = getc(fp);
            assert(c == 'N');
            c = getc(fp);
            assert(c == 'A');
            c = getc(fp);
            assert(c == 'R');
            c = getc(fp);
            assert(c == 'Y');
            binary = true;
        }
    }
    assert(c != EOF);
    int n_read;
    if (binary) {
        read_unsigned(&n_vertices, fp);
        read_unsigned(&n_faces, fp);
        read_unsigned(&n_edges, fp);
    } else {
        char line[1000];
        do
        {
            REQUIRE(fgets(line, 1000, fp) != 0);
        } while ((line[0] == '#') || (is_empty_line(line)));
        n_read = sscanf(line, "%u %u %u", &n_vertices, &n_faces,
                        &n_edges);
        REQUIRE(n_read == 3);
    }
}


template <class T>
void tIndexedTriangleSet<T>::read_off_vertex_ascii(FILE *fp, bool has_colors)
{
    value_type x, y, z;
    int n_read;
    n_read = treal<T>::scan(fp, &x);
    REQUIRE(n_read == 1);
    n_read = treal<T>::scan(fp, &y);
    REQUIRE(n_read == 1);
    n_read = treal<T>::scan(fp, &z);
    REQUIRE(n_read == 1);
    _vertices.push_back(Point3 (x, y, z));
    if (has_colors) {
        value_type r, g, b, a;
        n_read = treal<T>::scan(fp, &r);
        REQUIRE(n_read == 1);
        n_read = treal<T>::scan(fp, &g);
        REQUIRE(n_read == 1);
        n_read = treal<T>::scan(fp, &b);
        REQUIRE(n_read == 1);
        n_read = treal<T>::scan(fp, &a);
        REQUIRE(n_read == 1);
        _vertex_colors.push_back(ColorRgb (r, g, b));
    }
}


template <class T>
void tIndexedTriangleSet<T>::read_off_vertex_binary(FILE *fp, bool has_colors)
{
    Point3 v;
    v.read(fp);
    _vertices.push_back(v);
    if (has_colors) {
        ColorRgb c;
        c.read(fp);
        _vertex_colors.push_back(c);
        float a;
        read_float(&a, fp);
    }
}


template <class T>
void tIndexedTriangleSet<T>::read_off_face_ascii(
    FILE *fp,
    bool subdivision_enabled,
    value_type max_area)
{
    unsigned n_sides;
    int n_read;
    n_read = fscanf(fp, "%u", &n_sides);
    REQUIRE(n_read == 1);
    assert(n_sides >= 3);
    vector<unsigned> sides(n_sides);
    for (unsigned j = 0; j < n_sides; j++) {
        n_read = fscanf(fp, "%u", &sides[j]);
        REQUIRE(n_read == 1);
//        --sides[j]; // Hack to read ""bad"" files with index that begins at 1
        assert(sides[j] < _vertices.size());
    }
    // skip to the end of line
    {
        char line[1000];
        fgets(line, 1000, fp);
    }

    insert_face(sides, subdivision_enabled, max_area);
}


template <class T>
void tIndexedTriangleSet<T>::read_off_face_binary(
    FILE *fp,
    bool subdivision_enabled,
    value_type max_area)
{
    unsigned n_sides;
    read_unsigned(&n_sides, fp);
    assert(n_sides >= 3);
    vector<unsigned> sides(n_sides);
    for (unsigned j = 0; j < n_sides; j++) {
        read_unsigned(&sides[j], fp);
        assert(sides[j] < _vertices.size());
    }

    // @@@ only binary has this?
    int color_components = 0;
    read_int(&color_components, fp);
    assert(color_components == 0);

    insert_face(sides, subdivision_enabled, max_area);
}


template <class T>
void tIndexedTriangleSet<T>::read_off(
    FILE *fp,
    bool subdivision_enabled,
    value_type max_area)
{
    assert(fp);
    assert(!subdivision_enabled || (treal<T>::is_positive(max_area)));

    // read header
    bool has_colors = false;
    unsigned n_vertices = 0;
    unsigned n_faces = 0;
    unsigned n_edges = 0;
    bool binary = false;
    read_off_header(fp, has_colors, n_vertices, n_faces, n_edges, binary);

    if (binary) {
        // read vertices
        for (unsigned i = 0; i < n_vertices; i++) {
            read_off_vertex_binary(fp, has_colors);
        }
		
        // read faces
        for (unsigned i = 0; i < n_faces; i++) {
            read_off_face_binary(fp, subdivision_enabled,
                                 max_area);
        }
    } else {
        // read vertices
        for (unsigned i = 0; i < n_vertices; i++) {
            read_off_vertex_ascii(fp, has_colors);
        }

        // read faces
        for (unsigned i = 0; i < n_faces; i++) {
            read_off_face_ascii(fp, subdivision_enabled, max_area);
        }
    }

    if ((num_triangles() > 0) && (!has_vertex_normals())) {
        compute_vertex_normals();
        assert(has_vertex_normals());
    }
}


template <class T>
void tIndexedTriangleSet<T>::write_off_ascii(FILE *fp) const
{
    assert(fp);
    if (_vertex_colors.size() == 0) {
        fprintf(fp, "OFF\n");
    } else {
        fprintf(fp, "COFF\n");
    }
    fprintf(fp, "%d %d 0\n", _vertices.size(), _triangles.size());
    for (unsigned i = 0; i < _vertices.size(); i++) {
        const Point3 &p = _vertices[i];
        fprintf(fp, "%f %f %f", p.x(), p.y(), p.z());
        if (_vertex_colors.size() > 0) {
            const ColorRgb &c = _vertex_colors[i];
            fprintf(fp, " %f %f %f 1", c.r(), c.g(), c.b());
        }
        fprintf(fp, "\n");
    }
    for (unsigned i = 0; i < _triangles.size(); i++) {
        fprintf(fp, "3 %d %d %d\n", 
                _triangles[i].A(),
                _triangles[i].B(),
                _triangles[i].C());
    }
}


template <class T>
void tIndexedTriangleSet<T>::read(const char *file_name)
{
    FILE *fp = xfopen(file_name, "rb");
    char extension[10];
    get_file_extension(file_name, extension, sizeof(extension));
    if (!strcmp(extension, "off")) {
        read_off(fp);
    } else if (!strcmp(extension, "rtpi")) {
        read_rtpi(fp);
    } else {
        fprintf(stderr, "%s: unsupported model type\n", file_name);
        exit(EXIT_FAILURE);
    }
    fclose(fp);

    if ((num_triangles() > 0) && (!has_vertex_normals())) {
        compute_vertex_normals();
    }
}


template <class T>
bool tIndexedTriangleSet<T>::is_valid() const
{
    for (unsigned i = 0; i < num_triangles(); i++) {
        const IndexedTriangle &t = indexed_triangle(i);
        if ((t.A() < 0) || (t.A() >= (int) num_vertices()) ||
            (t.B() < 0) || (t.B() >= (int) num_vertices()) ||
            (t.C() < 0) || (t.C() >= (int) num_vertices())) {
            return false;
        }
    }
    return true;
}


template <class T>
void tIndexedTriangleSet<T>::add(const tIndexedTriangleSet &model)
{
    assert(model.is_valid());

    // remember old number of vertices
    unsigned n = num_vertices();

    // add vertices
    for (unsigned i = 0; i < model.num_vertices(); i++) {
        _vertices.push_back(model.vertex(i));
    }

    // add colors
    for (unsigned i = 0; i < model.num_vertex_colors(); i++) {
        _vertex_colors.push_back(model.vertex_color(i));
    }

    // add normals
    _vertex_normals.insert(
        _vertex_normals.end(),
        model._vertex_normals.begin(),
        model._vertex_normals.end());

    // add triangles
    for (unsigned i = 0; i < model.num_triangles(); i++) {
        const IndexedTriangle &t = model.indexed_triangle(i);
        _triangles.push_back(IndexedTriangle(t.A() + n,
                                             t.B() + n,
                                             t.C() + n));
    }

    assert(is_valid());
}


template <class T>
void tIndexedTriangleSet<T>::transform(const Matrix4 &m)
{
    for (unsigned i = 0; i < _vertices.size(); i++) {
        _vertices[i].transform(m);
    }

    // See Foley et al., p. 217
    Matrix4 m2 = m.inverse().transpose();

    for (unsigned i = 0; i < _vertex_normals.size(); i++) {
        _vertex_normals[i].transform(m2);
    }
}


template <class T>
void tIndexedTriangleSet<T>::map_colors(const Image &img, const Camera &cam)
{
    if (_vertex_colors.size() == 0) {
        _vertex_colors.resize(_vertices.size());
    }
    Viewport vp(0, 0, img.width(), img.height());
    tView<T> view(cam, vp);
    for (unsigned i = 0; i < _vertices.size(); i++) {
        const Point3 &p = _vertices[i];
        if (cam.sees(p)) {
            Point2 q = view.viewport_point(p);
            int qx = (int) (q.x() + 0.5);
            int qy = (int) (q.y() + 0.5);
            if ((qx >= 0) && (qx < img.width())
                && (qy >= 0) && (qy < img.height())) {
                // In theory we don't need the above
                // if, but in practice we do, because
                // of rounding errors.
                _vertex_colors[i] = img.pixel_rgb(qx, qy);
            }
        } else {
            _vertex_colors[i] = COLOR_RGB_GRAY50;
        }
    }
}


template <class T>
void tIndexedTriangleSet<T>::add_vertex(const Point3 &p)
{
    _vertices.push_back(p);
}


template <class T>
void tIndexedTriangleSet<T>::add_vertex(const Point3 &p, const ColorRgb &c)
{
    _vertices.push_back(p);
    _vertex_colors.push_back(c);
}


template <class T>
void tIndexedTriangleSet<T>::subdivide(const IndexedTriangle &it, value_type max_area)
{
    bool has_colors = _vertex_colors.size() > 0;
    bool has_normals = _vertex_normals.size() > 0;

    Point3 pa = _vertices[it.A()];
    Point3 pb = _vertices[it.B()];
    Point3 pc = _vertices[it.C()];
    Triangle3 t(pa, pb, pc);
    if (t.area() > max_area) {
        Point3 ab = Point3::midpoint(t.A(), t.B());
        Point3 ac = Point3::midpoint(t.A(), t.C());
        Point3 bc = Point3::midpoint(t.B(), t.C());

        int iab = _vertices.size();
        _vertices.push_back(ab);
        int iac = _vertices.size();
        _vertices.push_back(ac);
        int ibc = _vertices.size();
        _vertices.push_back(bc);

        IndexedTriangle t1(it.A(), iab, iac);
        IndexedTriangle t2(it.B(), ibc, iab);
        IndexedTriangle t3(it.C(), iac, ibc);
        IndexedTriangle t4(iab, ibc, iac);

        if (has_colors) {
            _vertex_colors.push_back(
                0.5 * (_vertex_colors[it.A()] +
                       _vertex_colors[it.B()]));
            _vertex_colors.push_back(
                0.5 * (_vertex_colors[it.A()] +
                       _vertex_colors[it.C()]));
            _vertex_colors.push_back(
                0.5 * (_vertex_colors[it.B()] +
                       _vertex_colors[it.C()]));
        }

        if (has_normals) {
            Vector3 npa = (_vertex_normals[it.A()] + 
                           _vertex_normals[it.B()]) * (T)0.5;
            Vector3 npb = (_vertex_normals[it.A()] + 
                           _vertex_normals[it.C()]) * (T)0.5;
            Vector3 npc = (_vertex_normals[it.B()] + 
                           _vertex_normals[it.C()]) * (T)0.5;
            _vertex_normals.push_back(npa);
            _vertex_normals.push_back(npb);
            _vertex_normals.push_back(npc);
        }
        subdivide(t1, max_area);
        subdivide(t2, max_area);
        subdivide(t3, max_area);
        subdivide(t4, max_area);
    } else {
        _triangles.push_back(it);
        if (has_normals) {		
            _vertex_normals.push_back(t.normal());
        }
    }

    assert((_vertex_colors.size() == 0) ||
           (_vertex_colors.size() == _vertices.size()));
    assert((_vertex_normals.size() == 0) ||
           (_vertex_normals.size() == _vertices.size()));
}


typedef std::vector<unsigned> adjacent_face_list_t;


class vertex_adjacency_list {
public:
    explicit vertex_adjacency_list(unsigned n_vertices);
    ~vertex_adjacency_list();
    adjacent_face_list_t &operator[](unsigned i);

private:
    vector<adjacent_face_list_t *> _adj_list;
};


vertex_adjacency_list::vertex_adjacency_list(unsigned n_vertices)
    : _adj_list(n_vertices)
{
    assert(_adj_list.size() == n_vertices);
    for (unsigned i = 0; i < n_vertices; i++) {
        _adj_list[i] = new adjacent_face_list_t;
    }
}


vertex_adjacency_list::~vertex_adjacency_list()
{
    for (unsigned i = 0; i < _adj_list.size(); i++) {
        delete _adj_list[i];
        _adj_list[i] = NULL;
    }
}


inline adjacent_face_list_t &vertex_adjacency_list::operator[](unsigned i)
{
    assert(i < _adj_list.size());
    assert(_adj_list[i] != NULL);
    return *(_adj_list[i]);
}


#if 0 // OLD FOR REFERENCE ONLY
// Compute the normal of vertices: vertex_normal = average of the
// normal to the adjacent faces.
template <class T>
void tIndexedTriangleSet<T>::compute_vertex_normals()
{
    vertex_adjacency_list adj_list(num_vertices());

    _vertex_normals.clear();
    _vertex_normals.reserve(num_vertices());

    // Build adjacency list.  For every vertex, remember it's
    // adjacent faces.
    for (unsigned f_index = 0; f_index < num_triangles(); ++f_index) {
        adj_list[_triangles[f_index].A()].push_back(f_index);
        adj_list[_triangles[f_index].B()].push_back(f_index);
        adj_list[_triangles[f_index].C()].push_back(f_index);
    }

    // Compute a normal per vertex.
    for (unsigned v_index = 0; v_index < num_vertices(); ++v_index) {
        adjacent_face_list_t &adj_faces = adj_list[v_index];

        unsigned n_adj_faces = adj_faces.size();
        if (n_adj_faces == 0) {
            //GTB_WARNING("vertex with no adjacent face");
        }

        Vector3 n(0.0);
        for (unsigned j = 0; j < n_adj_faces; j++) {
            unsigned adj_face = adj_faces[j];
            n += (face_normal(adj_face) / n_adj_faces);
        }
        n.normalize();
        _vertex_normals.push_back(n);
    }
    assert(_vertex_normals.size() == _vertices.size());
}
#endif // 0

// Compute the normal of vertices: vertex_normal = average of the
// normal to the adjacent faces.
template <class T>
void tIndexedTriangleSet<T>::compute_vertex_normals()
{

    _vertex_normals.clear();
    _vertex_normals.resize(num_vertices(), Vector3(0,0,0));

    unsigned F = num_triangles();

    for (unsigned f = 0; f < F; ++f)
    {
        Vector3 normal = face_normal(f);

        const IndexedTriangle& tri = indexed_triangle(f);

        _vertex_normals[tri.A()] += normal;
        _vertex_normals[tri.B()] += normal;
        _vertex_normals[tri.C()] += normal;

        // Q&D hack to avoid problems with data files that contains
        // triangles such as (v1,v2,v3) and (v1,v3,v2) (i.e. back to back)
        if ((_vertex_normals[tri.A()].squared_length() < 1e-6)) _vertex_normals[tri.A()] = normal;
        if ((_vertex_normals[tri.B()].squared_length() < 1e-6)) _vertex_normals[tri.B()] = normal;
        if ((_vertex_normals[tri.C()].squared_length() < 1e-6)) _vertex_normals[tri.C()] = normal;
    }

    unsigned V = num_vertices();
    for (unsigned v = 0 ; v < V; ++v)
    {
        _vertex_normals[v].normalize();
//        assert(fabs(_vertex_normals[v].squared_length() - 1) < 0.001);
    }
}



template <class T>
typename  tIndexedTriangleSet<T>::Vector3 tIndexedTriangleSet<T>::face_normal(unsigned face_index) const
{
    assert(face_index < num_triangles());
    int a = indexed_triangle(face_index).A();
    int b = indexed_triangle(face_index).B();
    int c = indexed_triangle(face_index).C();
    const Point3 &A = vertex(a);
    const Point3 &B = vertex(b);
    const Point3 &C = vertex(c);
    Vector3 normal = Point3::normal(A, B, C);

#if 0
    assert (isfinite(normal[0]));
    assert (isfinite(normal[1]));
    assert (isfinite(normal[2]));
#endif

    return normal;
}

template <class T>
void tIndexedTriangleSet<T>::scale(const Point3& origin, value_type a)
{
    typename vertex_list::iterator pf = _vertices.begin();
    typename vertex_list::iterator pl = _vertices.end();
    for (; pf != pl; ++pf)
    {
        *pf = pf->scale(origin, a);
    }
}

template <class T>
void tIndexedTriangleSet<T>::flip_normals()
{

    typename normal_list::iterator  pf = _vertex_normals.begin();
    typename normal_list::iterator  pl = _vertex_normals.end();
    for (; pf != pl; ++pf)
    {
        pf->flip();
    }

    for (unsigned f_index = 0; f_index < num_triangles(); ++f_index) {
        _triangles[f_index].flip();
    }
}

template <class T>
void tIndexedTriangleSet<T>::set_vertex(const Point3& p, unsigned idx)
{
    assert(idx<_vertices.size());

    _vertices[idx] = p;
}

template <class T>
void tIndexedTriangleSet<T>::set_vertex_color(unsigned idx, const ColorRgb& c)
{
    _vertex_colors[idx] = c;
}

template <class T>
void tIndexedTriangleSet<T>::need_colors()
{
    _vertex_colors.resize(_vertices.size());
}

template <class T>
void tIndexedTriangleSet<T>::need_face_colors()
{
    _face_colors.resize(_triangles.size());
}

template <class T>
bool tIndexedTriangleSet<T>::has_face_colors() const
{
    return _face_colors.size() > 0;
}

template <class T>
const typename tIndexedTriangleSet<T>::color_list &tIndexedTriangleSet<T>::face_colors() const
{
    return _face_colors;
}


/*
 * Compute the volume of the mesh, formula taken from 
 * "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow",
 *      Mathieu Desbrun Mark Meyer Peter Schr¨oder Alan H. Barr
 *
 *  V = 1/6 Sigma_{k = 1,n_faces} g_k * N_k
 *  g_k = (x^1_k + x^2_k + x^3_k)/3 - centroid of face
 *  N_x = normal to face
 */
template <class T>
typename  tIndexedTriangleSet<T>::value_type tIndexedTriangleSet<T>::volume() const
{
    value_type V = 0;
   
    unsigned F = num_triangles();

    for (unsigned k = 0; k < F; ++k)
    {
        const IndexedTriangle& tri = indexed_triangle(k);
        const Point3& p1 = vertex(tri.A());
        const Point3& p2 = vertex(tri.B());
        const Point3& p3 = vertex(tri.C());

//        Vector3 Nk = face_normal(k);
        Vector3 Nk = (p2-p1).cross(p3-p1);

        Point3 gk(0,0,0);


        gk.add(p1);
        gk.add(p2);
        gk.add(p3);

        gk.scalar_scale(1.0/3.0);


        V += Nk.dot(gk);
    }
    return V / 6.0;
}

template <class T>
void tIndexedTriangleSet<T>::remove_degenerated_faces()
{
    unsigned F = num_triangles();
    unsigned L = F;

    int NRemoved;
    do {
        NRemoved = 0;
        for (unsigned k = 0; k < F; ++k)
        {
            Vector3 n = face_normal(k);
            if (!(isfinite(n[0]) && isfinite(n[1]) && isfinite(n[2])))
            {
                printf("Removing face:%d\n", k);
                _triangles[k] = _triangles[L];
                --L;
                ++NRemoved;
                --k;
            }
        }
        _triangles.erase(_triangles.begin() + L, _triangles.end());
    } while (NRemoved);
}

/*
 * remove faces like <2 22 33> and <22 2 33>
 * removes both of the faces
 * report to stdout
 * method:
 *   sort faces and then look for the bad ones
 */
template <class T>
void tIndexedTriangleSet<T>::remove_double_sided_faces(bool remove_both)
{
    std::sort(_triangles.begin(), _triangles.end(), face_compare_less);

    unsigned F = num_triangles();

    std::vector<unsigned> eraselist;
    for (unsigned f = 0; f < F-1; ++f)
    {
        if (face_compare_eq(_triangles[f], _triangles[f+1]))
        {
            eraselist.push_back(f);
            if (remove_both) eraselist.push_back(f+1);
        }
    }

    if (eraselist.size()) printf("Removing duplicated triangles: ");
    std::vector<unsigned>::reverse_iterator f = eraselist.rbegin();
    std::vector<unsigned>::reverse_iterator l = eraselist.rend();
    for (; f != l; ++f)
    {
        printf(" %d", *f);
        _triangles.erase(_triangles.begin() + *f);
    }
    if (eraselist.size()) printf("\n");
}

/* 
 * A helper function to the above
 * return true of f1 < f2, false otherwise
 */
template <class T>
bool tIndexedTriangleSet<T>::face_compare_less(const IndexedTriangle& ct1, const IndexedTriangle& ct2)
{
    IndexedTriangle t1 = ct1;
    IndexedTriangle t2 = ct2;

    std::sort(t1.get_indices(), t1.get_indices()+3);
    std::sort(t2.get_indices(), t2.get_indices()+3);
    if (t1[0] == t2[0])
    {
        if (t1[1] == t2[1])
        {
            return t1[2] < t2[2];
        }
        else return t1[1] < t2[1];
    }
    else return t1[0] < t2[0];
}

template <class T>
bool tIndexedTriangleSet<T>::face_compare_eq(const IndexedTriangle& ct1, const IndexedTriangle& ct2)
{
    IndexedTriangle t1 = ct1;
    IndexedTriangle t2 = ct2;

    std::sort(t1.get_indices(), t1.get_indices()+3);
    std::sort(t2.get_indices(), t2.get_indices()+3);
    if ((t1[0] == t2[0]) && (t1[1] == t2[1]) && (t1[2] == t2[2])) return true;
    else return false;
}

/*
 * Remove vertices that do not belong to any triangle
 */
template <class T>
void tIndexedTriangleSet<T>::remove_lonely_vertices()
{
    unsigned V = num_vertices();
    std::vector<bool> lonelylist(V, true);

    unsigned F = num_triangles();
    for (unsigned f = 0; f < F; ++f)
    {
        const IndexedTriangle& tri = _triangles[f];
        for (unsigned j = 0; j < 3; ++j)
            lonelylist[tri[j]] = false;
    }

    for (int i = V-1; i >= 0; --i)
    {
        if (lonelylist[i])
        {
            _vertices.erase(_vertices.begin() + i);
            erase_vertex(i);
            printf("Erasing vertex: %d\n", i);
        }
    }
}

/*
 * Erase a vertex and update the face list
 * Assume the vertex does not belong to any face
 * Terrible O(N) algorithm
 */
template <class T>
void tIndexedTriangleSet<T>::erase_vertex(unsigned v)
{
    unsigned F = num_triangles();

    for (unsigned f = 0; f < F; ++f)
    {
        IndexedTriangle& tri = _triangles[f];
        for (int j = 0; j < 3; ++j)
        {
            if ((unsigned)tri[j] > v) --tri[j];
        }
    }
}

/*
 * Merge two vertices v1,v2 to v1
 * this means that any face that points to v1 will point to v2
 * This routine also, erases v2 properly, a bit slow, but working...
 */
template <class T>
void tIndexedTriangleSet<T>::edge_collaps(unsigned v1, unsigned v2)
{
    unsigned F = num_triangles();

    std::vector<unsigned> faceeraselist;

    for (unsigned f = 0; f < F; ++f)
    {
        IndexedTriangle& tri = _triangles[f];
        int edge_in_triangle = 0;
        for (int j = 0; j < 3; ++j)
        {
            if ((unsigned)tri[j] == v1)
            {
                ++edge_in_triangle;
            }
            if ((unsigned)tri[j] == v2) 
            {
//                tri[j] = v1;
                ++edge_in_triangle;
            }
        }
        if (edge_in_triangle == 2)
        {
            faceeraselist.push_back(f);
        }
    }

    if (faceeraselist.size() == 0) return;

    for (unsigned i = 0; i < faceeraselist.size(); ++i)
    {
        _triangles.erase(_triangles.begin() + faceeraselist[i]);        
    }

    _vertices.erase(_vertices.begin() + v2);

    for (unsigned f = 0; f < F; ++f)
    {
        IndexedTriangle& tri = _triangles[f];
        int edge_in_triangle = 0;
        for (int j = 0; j < 3; ++j)
        {
            if ((unsigned)tri[j] == v2) 
            {
                tri[j] = v1;
                ++edge_in_triangle;
            }
            else if ((unsigned)tri[j] > v2) --tri[j];
        }
    }
}

template <class T>
void tIndexedTriangleSet<T>::erase_face(unsigned f)
{
    _triangles.erase(_triangles.begin() + f);
}

template <class T>
void tIndexedTriangleSet<T>::apply(const Matrix4& M)
{
    typename vertex_list::iterator f = _vertices.begin();
    typename vertex_list::iterator l = _vertices.end();

    for (; f != l; ++f)
    {
        *f = M*(*f);
    }

    tModel<T>::invalidate_all();
}

template <class T>
void tIndexedTriangleSet<T>::clear()
{
    _vertices.clear();
    _vertex_normals.clear();
    _vertex_colors.clear();
    _face_colors.clear();
    _triangles.clear();

    tModel<T>::invalidate_all();
}

template class tIndexedTriangleSet<float>;
template class tIndexedTriangleSet<double>;

GTB_END_NAMESPACE
