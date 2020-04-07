
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


#ifndef __GTB_MEMORY_H
#define __GTB_MEMORY_H

#include "ptrs.h"
#include "counter.h"
#include "timer.h"
#include "hirestimer.h"

#ifdef WIN32
/*
 * Overload the default global new / delete operators
 * to use aligned allocation
 */
void * operator new(size_t size);
void operator delete(void * mem);
void * operator new[](size_t size);
void operator delete[](void * mem);

void start_leak_collection();
void stop_leak_collection();
void dump_leaks();
#endif

#endif // __GTB_MEMORY_H
