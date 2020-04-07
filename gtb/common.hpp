
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


#ifndef GTB_NAMESPACE_INCLUDED
#define GTB_NAMESPACE_INCLUDED

#define GTB_BEGIN_NAMESPACE namespace gtb {
#define GTB_END_NAMESPACE }

#ifdef WIN32
#include <gtb/windows.hpp>
#else
#include <gtb/unix.hpp>
#endif

#ifdef WIN32
#define BREAK __asm int 3
#else
#include <signal.h>
#define BREAK kill(getpid(), SIGQUIT)
#endif

#if defined(REAL_IS_FLOAT)
  #define GTB_GENERATE_STD_TYPEDEF(c) typedef c##f c
#else
  #define GTB_GENERATE_STD_TYPEDEF(c) typedef c##d c
#endif

#define GTB_GENERATE_CLASS_TYPEDEFS(c)		\
typedef t##c<float> c##f;			\
typedef t##c<double> c##d;			\
GTB_GENERATE_STD_TYPEDEF(c);

#endif // GTB_NAMESPACE_INCLUDED
