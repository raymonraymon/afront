
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


#ifndef GTB_SPHERE_INCLUDED
#define GTB_SPHERE_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/point3.hpp>
#include <gtb/graphics/vector3.hpp>
#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

GTB_BEGIN_NAMESPACE

template <class T>
class tSphere {
public:
    typedef T value_type;
    typedef tPoint3<T> Point3;

    tSphere();
    tSphere(const tSphere &s);
    tSphere(const Point3 &center, value_type radius);
    tSphere &operator=(const tSphere &s);

    bool operator==(const tSphere &s) const;
    bool operator!=(const tSphere &s) const;
    
    const Point3 &center() const;
    Point3 get_center() const;
    void set_center(const Point3& center);

    value_type radius() const;
    value_type get_radius() const;
    void set_radius(value_type r);

    bool contains(const Point3 &p) const;
    bool contains_approx(const Point3 &p) const;

    void render(int slices = 15,
                int stacks = 10,
                GLenum draw_style = GLU_LINE,
                GLenum normal_type = GLU_SMOOTH) const;

protected:
    Point3 _center;
    value_type _radius;
};

typedef tSphere<float> Spheref;
typedef tSphere<double> Sphered;
#if defined(REAL_IS_FLOAT)
typedef Spheref Sphere;
#else
typedef Sphered Sphere;
#endif


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/sphere.ipp>
#endif

#endif // GTB_SPHERE_INCLUDED
