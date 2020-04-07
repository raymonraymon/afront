
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


// Warnings

#pragma warning( disable : 4267 )
#pragma warning( disable : 4996 )
#pragma warning( disable : 4530 4652 4652 4651) // Exception handling

// OpenGL stuff

//#define WIN32_LEAN_AND_MEAN

#include <windows.h>

#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glext.h>

#ifndef GL_BGR
#define GL_BGR GL_BGR_EXT
#endif
#ifndef GL_BGRA
#define GL_BGRA GL_BGRA_EXT
#endif

// C stuff

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <io.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <malloc.h>

#ifdef min
#undef min
#undef max
#endif

#define S_ISREG(mode) (mode & _S_IFREG)
#define S_ISDIR(mode) (mode & _S_IFDIR)

#ifndef PATH_MAX
#define PATH_MAX 512
#endif

typedef int socklen_t;

// C++ stuff

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <stack>
#include <queue>

#define for if (0); else for

#define REAL_IS_DOUBLE 1

#define hypot _hypot
#define copysign _copysign

#define sleep Sleep

#define HAVE_STRICMP 1


