
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
#include "PlaneTransformation.h"

GTB_BEGIN_NAMESPACE

template <class REAL>
void plane_transformation<REAL>::LoadFromPlaneOGLMatrix() const
{
    GLdouble m[16];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if ((i==3) || (j==3)) m[i+j*4] = 0;
            else m[i+j*4] = Rt[i][j];
        
    m[15] = 1;
    glMultMatrixd(m);
    glTranslated(-T[0], -T[1], -T[2]);
}

template <>
void plane_transformation<float>::write(FILE* f) const
{
    m_plane.write(f);
    write_float(r[0], f);
    write_float(r[1], f);
    write_float(r[2], f);
    write_bool(_have_right_vector, f);
}

template <>
void plane_transformation<double>::write(FILE* f) const
{
    m_plane.write(f);
    write_double(r[0], f);
    write_double(r[1], f);
    write_double(r[2], f);
    write_bool(_have_right_vector, f);
}
    
template <>
void plane_transformation<float>::read(FILE* f)
{
    m_plane.read(f);
    read_float(&r[0], f);
    read_float(&r[1], f);
    read_float(&r[2], f);
    read_bool(&_have_right_vector, f);
    ComputeTransformation();
}

template <>
void plane_transformation<double>::read(FILE* f)
{
    m_plane.read(f);
    read_double(&r[0], f);
    read_double(&r[1], f);
    read_double(&r[2], f);
    read_bool(&_have_right_vector, f);
    ComputeTransformation();
}

template plane_transformation<float>;
template plane_transformation<double>;

GTB_END_NAMESPACE
