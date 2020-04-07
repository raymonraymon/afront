
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


#ifndef __GTB_PLANE_TRANSFORMATION_HPP
#define __GTB_PLANE_TRANSFORMATION_HPP

#include <gtb/common.hpp>
#include <gtb/graphics/plane.hpp>
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/surfelset.hpp>
#include <gtb/graphics/mat3.hpp>
#include <gtb/graphics/ogltools.h>

GTB_BEGIN_NAMESPACE


/*
 * Name: PlaneTransformation.hpp
 * Author: Shachar Fleishman
 * Description:
 *
 * A class that help transform points in 3D world coordinates to
 * and from 2D points in the coordinates space on a plane.
 *
 *   
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */

/*
 */
template<class REAL>
class plane_transformation
{
public:
    plane_transformation() : _have_right_vector(false) {}

    plane_transformation(const tPlane<REAL>& plane, const tPoint3<REAL>& center) :
        m_plane(plane), 
        r(center), 
        _have_right_vector(false)
    {
        ComputeTransformation();
    }

    void Set(const tPlane<REAL>& plane, const tPoint3<REAL>& center)
    {
        m_plane = plane;
        r = center;
        _have_right_vector = false;
        ComputeTransformation();
    }

    //
    // Set the value of the vector pointing right (orthogonal
    // to the normal to the plane) and recompute the transfromation
    // Assume:
    //    right is normalized.
    //    dot(right, normal) = 0
    //
    // The coordinate system is: n(normal), right, f(orward)
    //
    void SetRight(const tVector3<REAL>& right)
    {
        _right = right;
        _have_right_vector = true;
        ComputeTransformation();
    }

    //
    // Flip the orientation of the plane
    // i.e. make the positive and negative half-spaces flip
    // position
    //
    void FlipOrientation(bool value=true)
    {
        tVector3<REAL> right;
        if (_have_right_vector)
        {
            right = _right;
        }
        else
        {
            right = PlaneNormal2World(tVector3<REAL>(1,0,0));
        }
        right.flip();
        SetRight(right);
        m_plane.flip_orientation();

        ComputeTransformation();// NEED faster procedure
    }

    //
    //
    //
    tPoint3<REAL> FromPlane(const Point2& x) const
    {
        tPoint3<REAL> x3(x[0], x[1], 0);
        return Rt * (x3-T);
    }

    //
    //
    //
    tPoint3<REAL> FromPlane(const tPoint3<REAL>& x) const
    {
        return Rt * (x-T);
    }

    //
    //
    //
    void LoadFromPlaneOGLMatrix() const;

    const tPlane<REAL>& get_plane() const { return m_plane; }

    //
    //
    //
    tPoint3<REAL> ToPlane(const tPoint3<REAL>& x) const
    {
        return R*x + T;
    }

    //
    // Input: 3d normal vector in plane coordinates
    // Return: 3d normal vector in world coordinates
    //
    tVector3<REAL> PlaneNormal2World(const tVector3<REAL>& n) const
    {
        return Rt * n;
    }

    //
    // Input: 3d normal vector in world coordinates
    // Return: 3d normal vector in plane coordinates
    //
    tVector3<REAL> WorldNormal2Plane(const tVector3<REAL>& n) const
    {
        return R * n;
    }

    const tmat3<REAL>& GetR() { return R; }
    const tVector3<REAL>& GetT() { return T; }
    const tPoint3<REAL>& get_origin() const { return r; }

    void write(FILE* f) const;
    void read(FILE* f);
    
    bool verify() const
    {
        return (isfinite(r[0]) && isfinite(m_plane.normal()[0]) && isfinite(m_plane.D()));
    }
private:
    tPlane<REAL> m_plane;
    tPoint3<REAL> r;            // Define the origin:
    // r+m_plane.t*m_plane.n is a point 
    // on the plane. that defines the origin.
    // Why this special definition? because
    // that's what we need for the MLS projection!
    tmat3<REAL> R;                // Rotation matrix
    tmat3<REAL> Rt;               // Inverse rotation matrix

    tVector3<REAL> T;           // Translation vector

    bool _have_right_vector;
    tVector3<REAL> _right;

    //
    // Check Householder transformation:
    //  Q = I - 2*n*n'/(n'*n)
    //  or simply if n is normalized Q=I-2*n*n'
    //  where n is a column vector
    //
    void ComputeTransformation()
    {
        if (_have_right_vector)
        {
            assert(fabs(_right*m_plane.normal()) < 1e-5);
            
            tVector3<REAL> f = m_plane.normal().cross(_right);
            R[0] = _right;
            R[1] = f;
            R[2] = m_plane.normal();
        }
        else
        {
            R = cs_change(m_plane.normal()).transpose();
            // Why (r + (p.n * p.D) and not simply r?
            // because the plane is relative to r, i.e. p.D is from r and not 0.
        }

        tVector3<REAL> vr = tVector3<REAL>(r);
        T = -(R * (vr + (m_plane.normal() * m_plane.D() )));

        Rt = R.transpose();
    }
};

typedef plane_transformation<float> plane_transformationf;
typedef plane_transformation<double> plane_transformationd;

//
// Transform a set of points in 3D to 2D
//
#if 0
template <class SURFELSET>
void surfel_set_to_2d_surfel_set (
    const SURFELSET& src,
    const plane_transformation& xform,
    gtb::point3_set& dest)
{
    typename SURFELSET::const_iterator f = src.begin();
    typename SURFELSET::const_iterator l = src.end();
    for (; f != l; ++f)
    {
        dest.push_back(xform.ToPlane(f->p));
    }
}

#endif

template <class VIEW, class T>
inline void view_to_plane(
    const VIEW& view,
    const plane_transformation<T>& xform,
    point3_set& dest)
{
    unsigned n_points = view.size();
    for (unsigned i = 0; i < n_points; ++i)
    {
        dest.push_back(xform.ToPlane(view.vertex(i)));
    }
}

GTB_END_NAMESPACE

#endif // __GTB_PLANE_TRANSFORMATION_HPP
