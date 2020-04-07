
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


using namespace std;

GTB_BEGIN_NAMESPACE

template<> const tPoint2<float> tPoint2<float>::ZERO(0.0f, 0.0f);
template<> const tPoint2<double> tPoint2<double>::ZERO(0.0, 0.0);

template <class T>
tPoint2<T> &tPoint2<T>::scale(const tPoint2 &origin, T s)
{
    tVector2<T> t = origin - ZERO;
    *this -= t;
    _p[0] *= s;
    _p[1] *= s;
    *this += t;
    return *this;
}


template <class T>
tPoint2<T> &tPoint2<T>::rotate(T theta)
{
    T c = cos(theta);
    T s = sin(theta);
    T px = _p[0];
    T py = _p[1];
    _p[0] = (px * c) - (py * s);
    _p[1] = (px * s) + (py * c);
    return *this;
}


template <class T>
tPoint2<T> &tPoint2<T>::rotate(const tPoint2 &origin, T theta)
{
    tVector2<T> t(origin._p);
    *this -= t;
    t.rotate(theta);
    *this += t;
    return *this;
}


template <class T>
tPoint2<T> tPoint2<T>::centroid(const vector<tPoint2> &v)
{
    T cx = 0.0;
    T cy = 0.0;
    for (unsigned i = 0; i < v.size(); i++) {
        cx += v[i].x() / v.size();
        cy += v[i].y() / v.size();
    }
    return tPoint2(cx, cy);
}


/*----------------------------------------------*/
template <class T>
T LeftTurn(const tPoint2<T>& p0, const tPoint2<T>& p1, const tPoint2<T>& p2) 
{
    return (p1.x()-p0.x())*(p2.y()-p1.y()) -
        (p2.x()-p1.x())*(p1.y()-p0.y());
}

template float LeftTurn(const tPoint2<float>& p0, const tPoint2<float>& p1, const tPoint2<float>& p2);
template double LeftTurn(const tPoint2<double>& p0, const tPoint2<double>& p1, const tPoint2<double>& p2);


/// by Olga - this function returns true if the given
/// point is inside a triangle defined by p0,p1,p2
/// Pre-condition: both point and *this are on xy plane
template <class T>
bool ContainsPoint(const tPoint2<T>& p0, const tPoint2<T>& p1, const tPoint2<T>& p2, const tPoint2<T>& point)
{
    T L0 = LeftTurn(p0, p1, point);
    T L1 = LeftTurn(p1, p2, point);
    T L2 = LeftTurn(p2, p0, point);

    return ((L0<=0 && L1<=0 && L2<=0) ||
            (L0>=0 && L1>=0 && L2>=0));
}

template <> bool ContainsPoint(const tPoint2<float>& p0, const tPoint2<float>& p1, const tPoint2<float>& p2, const tPoint2<float>& point);
template <> bool ContainsPoint(const tPoint2<double>& p0, const tPoint2<double>& p1, const tPoint2<double>& p2, const tPoint2<double>& point);

/*
 *  http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
 */
template <class T>
void circumcircle(const tPoint2<T>& p1, const tPoint2<T>& p2, const tPoint2<T>& p3, tPoint2<T>& center)
{
    T D = (p1.x() - p3.x())*(p2.y()-p3.y()) - (p2.x() - p3.x())*(p1.y() - p3.y());
    if (fabs(D) < 1e-8)
    {
        // This test is not debugged
        center = p1;
    }
    else
    {
        center[0] = 
            (((p1.x() - p3.x()) * (p1.x() + p3.x()) + (p1.y() - p3.y()) * (p1.y() + p3.y())) / 2 * (p2.y() - p3.y()) 
             -  ((p2.x() - p3.x()) * (p2.x() + p3.x()) + (p2.y() - p3.y()) * (p2.y() + p3.y())) / 2 * (p1.y() - p3.y())) 
            / D;

        center[1] = (((p2.x() - p3.x()) * (p2.x() + p3.x()) + (p2.y() - p3.y()) * (p2.y() + p3.y())) / 2 * (p1.x() - p3.x())
                     -  ((p1.x() - p3.x()) * (p1.x() + p3.x()) + (p1.y() - p3.y()) * (p1.y() + p3.y())) / 2 * (p2.x() - p3.x()))
            / D;
    }
}

template void circumcircle(const tPoint2<float>& p1, const tPoint2<float>& p2, const tPoint2<float>& p3, tPoint2<float>& center);
template void circumcircle(const tPoint2<double>& p1, const tPoint2<double>& p2, const tPoint2<double>& p3, tPoint2<double>& center);

template<>
void tPoint2<double>::read(FILE *fp)
{
	read_double(&_p[0], fp);
	read_double(&_p[1], fp);
}

template<>
 void tPoint2<double>::write(FILE *fp) const
{
	write_double(_p[0], fp);
	write_double(_p[1], fp);
}

template<>
void tPoint2<float>::read(FILE *fp)
{
	read_float(&_p[0], fp);
	read_float(&_p[1], fp);
}

template<>
void tPoint2<float>::write(FILE *fp) const
{
	write_float(_p[0], fp);
	write_float(_p[1], fp);
}


template class tPoint2<float>;
template class tPoint2<double>;

GTB_END_NAMESPACE
