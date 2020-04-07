
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
#include <gtb/graphics/line3.hpp>
#include <gtb/graphics/ray3.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/line3.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE

#if 0
const Line3 LINE3_POSITIVE_X(Point3(0.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0));
const Line3 LINE3_NEGATIVE_X(Point3(0.0, 0.0, 0.0), Vector3(-1.0, 0.0, 0.0));
const Line3 LINE3_POSITIVE_Y(Point3(0.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0));
const Line3 LINE3_NEGATIVE_Y(Point3(0.0, 0.0, 0.0), Vector3(0.0, -1.0, 0.0));
const Line3 LINE3_POSITIVE_Z(Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0));
const Line3 LINE3_NEGATIVE_Z(Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, -1.0));
#endif


/*-------------------------------------------------------*/
template<class T>
bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, T& sc, T& tc)
{
    tVector3<T> d12 = L1.origin() - L2.origin();
    T a = L1.direction().dot(d12);
    T b = L1.direction().squared_length();
    T c = L1.direction().dot(L2.direction());
    T d = L2.direction().dot(d12);
    // double e = c;
    T f = L2.direction().squared_length();
    T dt = f*b-c*c; // Determinant 
    if (fabs(dt) < 1e-8)
    {
        // Parallel llines
        sc = -d/c;
        tc = a/c;
    }
    else
    {
        // General lines
        tc = (b*d - a*c)/dt;
        sc = -(a - c*tc)/b;
    }
    return true;
}

template
bool LineDistance(const tLine3<float>& L1, const tLine3<float>& L2, float& sc, float& tc);
template
bool LineDistance(const tLine3<double>& L1, const tLine3<double>& L2, double& sc, double& tc);

template<class T>
bool LineDistance(const tLine3<T>& L1, const tLine3<T>& L2, tPoint3<T>& pc, tPoint3<T>& qc)
{
    T sc,tc;
    if (!LineDistance(L1, L2, sc, tc)) return false;

    pc = L1.point(sc);
    qc = L2.point(tc);

    return true;
}
template
bool LineDistance(const tLine3<float>& L1, const tLine3<float>& L2, tPoint3<float>& pc, tPoint3<float>& qc);
template
bool LineDistance(const tLine3<double>& L1, const tLine3<double>& L2, tPoint3<double>& pc, tPoint3<double>& qc);

#if 0
bool LineDistance(const tLine3& L1, const tLine3& L2, double& sc, double& tc)
{
    assert(fabs(L1._d.squared_length()-1.0) < 1e-8);
    assert(fabs(L2._d.squared_length()-1.0) < 1e-8);

    if (fabs(fabs(L1._d.dot(L2._d)) - 1) < 1e-8) return false;
    
    tVector3<T> w0 = L1._p-L2._p;
    // double a = L1._d.dot(L1._d);
    double b = L1._d.dot(L2._d);
    // double c = L2._d.dot(L2._d);
    double d = w0.dot(L1._d);
    double e = w0.dot(L2._d);
    sc = (b*e-/*c* */d)/(/* a*c */1-b*b);
    tc = (/*a* */e-b*d)/(/* a*c */1-b*b);

    return true;
}
#endif


template class tLine3<float>;
template class tLine3<double>;

GTB_END_NAMESPACE
