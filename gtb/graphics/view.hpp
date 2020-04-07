
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


#ifndef GTB_VIEW_INCLUDED
#define GTB_VIEW_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/viewport.hpp>
#include <gtb/graphics/ray3.hpp>


GTB_BEGIN_NAMESPACE


template <class T>
class tView {
public:
    typedef T value_type;
    typedef tCamera<T> Camera;
    typedef tVector3<T> Vector3;
    typedef tPoint3<T> Point3;
    typedef tPoint2<T> Point2;
    typedef tRay3<T> Ray3;
    typedef tBox3<T> Box3;

    tView();
    tView(const tView &view);
    tView(const Camera &cam, const Viewport &vp);
    tView &operator=(const tView &view);

    const Camera &camera() const;
    const Camera &get_camera() const;
    Camera &camera();
    Camera &get_camera();
    const Viewport &viewport() const;
    const Viewport &get_viewport() const;
    Viewport &viewport();
    Viewport &get_viewport();

    void set_camera(const Camera &cam);
    void set_viewport(const Viewport &vp);

    void rotate_world(const Point3 &origin, value_type dx, value_type dy);
    void translate_world(const Point3 &origin, value_type dx, value_type dy);
    void scale_world(const Point3 &origin, value_type dx, value_type dy);
    void scale_world(value_type factor);

    void load_camera() const;
    void gl_load_camera() const;
    void load_camera(const std::vector<Light> &headlights) const;
    void init_exterior_view(const Box3 &box);
    void init_interior_view(const Box3 &box);
    void init_image_view(unsigned width, unsigned height);

    Ray3 world_ray(int x, int y) const; // (0,0) is lower left corner
    Point3 viewport_point(int x, int y) const;
    Point2 viewport_point(const Point3 &p) const;

protected:
    Camera _camera;
    Viewport _viewport;
};

typedef tView<float> Viewf;
typedef tView<double> Viewd;
#if defined(REAL_IS_FLOAT)
typedef Viewf View;
#else
typedef Viewd View;
#endif


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/view.ipp>
#endif

#endif // GTB_VIEW_INCLUDED
