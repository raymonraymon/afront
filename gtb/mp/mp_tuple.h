
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


#ifndef __GTB_META_PROGRAMMING_TUPLE_H
#define __GTB_META_PROGRAMMING_TUPLE_H

#include <gtb/real/real.hpp>
#include <gtb/mp/mp.h>

GTB_BEGIN_NAMESPACE
GTB_MP_BEGIN_NAMESPACE

// These are more like lists than tuples, but anyway.

template <class T, int D> struct Tuple;

// Closest things to templated typedefs;
template <int D> struct IntTuple { typedef Tuple<int, D> Type; };
template <int D> struct RealTuple { typedef Tuple<real_t, D> Type; };

template <class T>
struct Tuple<T, 0> {
    Tuple() {};
};

template <class T, int D>
struct Tuple {
    Tuple(T v, Tuple<T, D-1> r): _value(v), _rest(r) {};
    T _value;
    Tuple<T, D-1> _rest;
};


// 9 dimensions should be enough for everyone
#define GTB_MP_MAKE_INT_TUPLE_1(a)                 (gtb::mp::Tuple<int, 1>((a), gtb::mp::Tuple<int,0>()))
#define GTB_MP_MAKE_INT_TUPLE_2(a,b)               (gtb::mp::Tuple<int, 2>((a), GTB_MP_MAKE_INT_TUPLE_1(b)))
#define GTB_MP_MAKE_INT_TUPLE_3(a,b,c)             (gtb::mp::Tuple<int, 3>((a), GTB_MP_MAKE_INT_TUPLE_2(b,c)))
#define GTB_MP_MAKE_INT_TUPLE_4(a,b,c,d)           (gtb::mp::Tuple<int, 4>((a), GTB_MP_MAKE_INT_TUPLE_3(b,c,d)))
#define GTB_MP_MAKE_INT_TUPLE_5(a,b,c,d,e)         (gtb::mp::Tuple<int, 5>((a), GTB_MP_MAKE_INT_TUPLE_4(b,c,d,e)))
#define GTB_MP_MAKE_INT_TUPLE_6(a,b,c,d,e,f)       (gtb::mp::Tuple<int, 6>((a), GTB_MP_MAKE_INT_TUPLE_5(b,c,d,e,f)))
#define GTB_MP_MAKE_INT_TUPLE_7(a,b,c,d,e,f,g)     (gtb::mp::Tuple<int, 7>((a), GTB_MP_MAKE_INT_TUPLE_6(b,c,d,e,f,g)))
#define GTB_MP_MAKE_INT_TUPLE_8(a,b,c,d,e,f,g,h)   (gtb::mp::Tuple<int, 8>((a), GTB_MP_MAKE_INT_TUPLE_7(b,c,d,e,f,g,h)))
#define GTB_MP_MAKE_INT_TUPLE_9(a,b,c,d,e,f,g,h,i) (gtb::mp::Tuple<int, 9>((a), GTB_MP_MAKE_INT_TUPLE_8(b,c,d,e,f,g,h,i)))
#define GTB_MP_MAKE_REAL_TUPLE_1(a)                 (gtb::mp::Tuple<gtb::real_t, 1>((a), gtb::mp::Tuple<gtb::real_t,0>()))
#define GTB_MP_MAKE_REAL_TUPLE_2(a,b)               (gtb::mp::Tuple<gtb::real_t, 2>((a), GTB_MP_MAKE_REAL_TUPLE_1(b)))
#define GTB_MP_MAKE_REAL_TUPLE_3(a,b,c)             (gtb::mp::Tuple<gtb::real_t, 3>((a), GTB_MP_MAKE_REAL_TUPLE_2(b,c)))
#define GTB_MP_MAKE_REAL_TUPLE_4(a,b,c,d)           (gtb::mp::Tuple<gtb::real_t, 4>((a), GTB_MP_MAKE_REAL_TUPLE_3(b,c,d)))
#define GTB_MP_MAKE_REAL_TUPLE_5(a,b,c,d,e)         (gtb::mp::Tuple<gtb::real_t, 5>((a), GTB_MP_MAKE_REAL_TUPLE_4(b,c,d,e)))
#define GTB_MP_MAKE_REAL_TUPLE_6(a,b,c,d,e,f)       (gtb::mp::Tuple<gtb::real_t, 6>((a), GTB_MP_MAKE_REAL_TUPLE_5(b,c,d,e,f)))
#define GTB_MP_MAKE_REAL_TUPLE_7(a,b,c,d,e,f,g)     (gtb::mp::Tuple<gtb::real_t, 7>((a), GTB_MP_MAKE_REAL_TUPLE_6(b,c,d,e,f,g)))
#define GTB_MP_MAKE_REAL_TUPLE_8(a,b,c,d,e,f,g,h)   (gtb::mp::Tuple<gtb::real_t, 8>((a), GTB_MP_MAKE_REAL_TUPLE_7(b,c,d,e,f,g,h)))
#define GTB_MP_MAKE_REAL_TUPLE_9(a,b,c,d,e,f,g,h,i) (gtb::mp::Tuple<gtb::real_t, 9>((a), GTB_MP_MAKE_REAL_TUPLE_8(b,c,d,e,f,g,h,i)))

GTB_MP_END_NAMESPACE
GTB_END_NAMESPACE

#endif
