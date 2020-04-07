
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
#include <gtb/math/math.hpp>
#endif // WIN32

GTB_BEGIN_NAMESPACE

#if 0
template <class T>
void spherical_to_euclidian(T phi, T theta, tVector3<T>& d)
{
    d[0] = (T)cos(phi)*cos(theta);
    d[1] = (T)sin(phi)*cos(theta);
    d[2] = (T)sin(theta);
}
template void spherical_to_euclidian<float>(float phi, float theta, tVector3<float>& d);
template void spherical_to_euclidian<double>(double phi, double theta, tVector3<double>& d);

template <class T>
void euclidian_to_spherical(const tVector3<T>& d, T& phi, T& theta)
{
    theta = (T)(asin(d[2]));
#if 0
    T cos_theta = cos(theta);
    T phi=asin(d[1]/cos_theta);
    T phi=acos(d[0]/cos_theta);
    phi = ac_phi;
    if (as_phi < 0)
    {
       phi = 2*M_PI+as_phi;
    }
    else if (ac_phi < 0) phi = ac_phi + M_PI;
#else
    phi = (T)(atan2(d[1], d[0]));
#endif
}
template void euclidian_to_spherical<float>(const tVector3<float>& d, float& phi, float& theta);
template void euclidian_to_spherical<double>(const tVector3<double>& d, double& phi, double& theta);
#endif

GTB_END_NAMESPACE

#ifdef OUTLINE
#define inline
#include <gtb/math/math.ipp>
#undef inline
#endif
