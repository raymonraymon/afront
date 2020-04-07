
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


#ifndef GTB_CAMERA_INCLUDED
#define GTB_CAMERA_INCLUDED

#include <gtb/common.hpp>
#include <gtb/math/math.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#include <gtb/graphics/coordinate_system.hpp>
#include <gtb/graphics/box3.hpp>
#include <gtb/graphics/plane.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/graphics/viewport.hpp>
#include <gtb/graphics/light.hpp>


GTB_BEGIN_NAMESPACE

template<class T>
class tCamera {
public:
    typedef T value_type;
    typedef tPoint3<T> Point3;
    typedef tVector3<T> Vector3;
    typedef tMatrix4<T> Matrix4;
    typedef tBox3<T> Box3;
    typedef tLine3<T> Line3;

    enum Type { PERSPECTIVE, ORTHOGRAPHIC };

    tCamera();

    tCamera(const Point3 &origin,
           const Vector3 &towards,
           const Vector3 &up,
           value_type x_fov,
           value_type y_fov,
           value_type near_distance,
           value_type far_distance);

    tCamera(const tCamera &cam);

    tCamera(const tCamera &cam,
           unsigned i,
           unsigned j,
           unsigned n_rows,
           unsigned n_cols);

    tCamera &operator=(const tCamera &cam);

    void reset(const Point3 &origin,
               const Vector3 &towards,
               const Vector3 &up);

    void reset(const Point3 &origin,
               const Vector3 &towards,
               const Vector3 &up,
               value_type x_fov,
               value_type y_fov,
               value_type near_distance,
               value_type far_distance);

    void reset_fov(value_type fovx, value_type fovy);

    const Point3 &origin() const;
    const Point3 &get_origin() const;
    Vector3 towards() const;
    const Vector3 &backwards() const;
    const Vector3 &up() const;
    const Vector3 &get_up() const;
    Vector3 down() const;
    const Vector3 &right() const;
    const Vector3 &get_right() const;
    Vector3 left() const;
    value_type x_fov() const;
    value_type y_fov() const;
    value_type aspect() const;
    value_type near_distance() const;
    value_type far_distance() const;
    const tCoordinateSystem<T> &coordinate_system() const;

    void load(const Viewport &viewport) const;
    void gl_load(const Viewport &viewport) const;
    void load(const Viewport &viewport,
              const std::vector<Light> &headlights) const;
    //void render() const;
    void render(value_type scale) const;

//    void view_all(const Box3 &box);
    void init_exterior_view(const Box3 &box);
    void init_interior_view(const Box3 &box);
    void init_image_view();

    void rotate(const Line3 &l, value_type theta);
    void pitch(value_type angle);
    void yaw(value_type angle);
    void roll(value_type angle);

    void translate(const Vector3 &t);
    void move_to(const Point3 &p);

    bool sees(const Point3 &point) const;
    bool sees(const Box3 &box) const;
    bool sees_ignoring_near_plane(const Box3 &b) const;
    bool contains(const Box3 &box) const;

    bool read(FILE *fp);
    void write(FILE *fp) const;

    bool read(const char *file_name);
    void write(const char *file_name) const;

    Matrix4 matrix() const;
    Matrix4 inverse_matrix() const;

    Point3 world_to_camera(const Point3 &p) const;
    Point3 camera_to_world(const Point3 &p) const;

    template <typename T2>
    friend std::ostream &operator<<(std::ostream &s, const tCamera<T2> &c);

protected:

    typedef tPlane<T> Plane;

    void setup_viewport(const Viewport &viewport) const;
    void setup_frustum() const;
    void load_headlights(const std::vector<Light> &headlights) const;
    void setup_look_at() const;

    tCoordinateSystem<T> _cs;
    value_type _x_fov, _y_fov;
    value_type _near_distance, _far_distance;
    Plane _planes[6];
    Point3 _corners[8];
    unsigned _row, _col;
    unsigned _n_rows, _n_cols;
    Type _type;
};

typedef tCamera<float> Cameraf;
typedef tCamera<double> Camerad;
#if defined(REAL_IS_FLOAT)
typedef Cameraf Camera;
#else
typedef Camerad Camera;
#endif

template <class T>
std::ostream &operator<<(std::ostream &s, const tCamera<T> &c);


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/camera.ipp>
#endif

#endif // GTB_CAMERA_INCLUDED
