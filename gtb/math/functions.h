
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


#ifndef __GTB_FUNCTIONS_H
#define __GTB_FUNCTIONS_H

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

//
// Some functions that I like
//

//
// A function that is
//  f(0) = 1
//  f(1) = 0
//  f'(0) = f'(1) = 0
template <class T> T bounded_decay(T x);

//
// expands bounded_decay to be defined in [-inf, inf]
//
template <class T> T infinit_bounded_decay(T x);

//
// Generalized inifite bounded decay
// so that the range and position of the decay is set by the user
//
// t - "time" of center of decay - i.e. how much to shift the function 
// s - scale of the decay, that is the decay is in the range [t-s, t+s]
template <class T> T infinit_bounded_decay(T x, T t, T s);

//
// inverted infinit_bounded_decay 
// that is equal to 0 at 0, 1 at 1 (or scaled)
//
template <class T> T inverse_infinit_bounded_decay(T x);
template <class T> T inverse_infinit_bounded_decay(T x, T t, T s);

/*------------------ Implementations -------------*/
template <class T> inline T bounded_decay(T x)
{
    assert( (x>=0) && (x<=1));
    return 2*x*x*x - 3*x*x + 1;
}

template <class T> inline T infinit_bounded_decay(T x)
{
    if (x < 0) return 1;
    else if (x > 1) return 0;
    else return bounded_decay(x);
}

template <class T> inline T infinit_bounded_decay(T x, T t, T s)
{
    return infinit_bounded_decay( (x-t)/(2*s) + 0.5);
}

template <class T> inline T inverse_infinit_bounded_decay(T x)
{
    return 1 - infinit_bounded_decay(x);
}

template <class T> inline T inverse_infinit_bounded_decay(T x, T t, T s)
{
    return 1 - infinit_bounded_decay(x, t, s);
}


GTB_END_NAMESPACE

#endif // __GTB_FUNCTIONS_H
