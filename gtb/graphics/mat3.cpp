
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
#include <gtb/graphics/point2.hpp>
#include <gtb/graphics/vector2.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/point2.ipp>
#undef inline
#endif

GTB_BEGIN_NAMESPACE

//mat3 mat3I(Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1));
//mat3 mat3Zero(0);

template <class T>
void tmat3<T>::FromEulerAnglesXYZ(value_type fYAngle, value_type fPAngle, value_type fRAngle)
{
    value_type fCos, fSin;

    fCos = cos(fYAngle);
    fSin = sin(fYAngle);
    tmat3 kXMat(
        Vector3(1.0,0.0,0.0),
        Vector3(0.0,fCos,-fSin),
        Vector3(0.0,fSin,fCos));

    fCos = cos(fPAngle);
    fSin = sin(fPAngle);
    tmat3 kYMat(
        Vector3(fCos,0.0,fSin),
        Vector3(0.0,1.0,0.0),
        Vector3(-fSin,0.0,fCos));

    fCos = cos(fRAngle);
    fSin = sin(fRAngle);
    tmat3 kZMat(
        Vector3(fCos,-fSin,0.0),
        Vector3(fSin,fCos,0.0),
        Vector3(0.0,0.0,1.0));

    *this = kXMat*(kYMat*kZMat);    
}

/*
 * Return:
 *   true if the solutin is unique
 */
template <class T>
bool tmat3<T>::ToEulerAnglesXYZ(value_type& rfXAngle, value_type& rfYAngle, value_type& rfZAngle) const
{
    // rot =  cy*cz          -cy*sz           sy
    //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
    //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

    if ( v[0][2] < (value_type)1.0 )
    {
        if ( v[0][2] > -(value_type)1.0 )
        {
            rfXAngle = atan2(-v[1][2],v[2][2]);
            rfYAngle = (value_type)asin((double)v[0][2]);
            rfZAngle = atan2(-v[0][1],v[0][0]);
            return true;
        }
        else
        {
            // WARNING.  Not unique.  XA - ZA = -atan2(r10,r11)
            rfXAngle = -atan2(v[1][0],v[1][1]);
            rfYAngle = M_PI_2;
            rfZAngle = (value_type)0.0;
            return false;
        }
    }
    else
    {
        // WARNING.  Not unique.  XAngle + ZAngle = atan2(r10,r11)
        rfXAngle = atan2(v[1][0],v[1][1]);
        rfYAngle = M_PI_2;
        rfZAngle = (value_type)0.0;
        return false;
    }    
}

template class tmat3<float>;
template class tmat3<double>;

GTB_END_NAMESPACE
