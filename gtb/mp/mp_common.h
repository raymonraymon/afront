
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


#ifndef __GTB_META_PROGRAMMING_COMMON_H
#define __GTB_META_PROGRAMMING_COMMON_H

#define GTB_MP_BEGIN_NAMESPACE namespace mp {
#define GTB_MP_END_NAMESPACE }

#ifdef WIN32
#define MP_INTEGER __int64
#else
#define MP_INTEGER long long int
#endif

#endif
