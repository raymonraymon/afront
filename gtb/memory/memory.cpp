
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
#include <gtb/memory/memory.hpp>
#endif // WIN32

#define USE_DL_PREFIX
#include "fmalloc.h"

#define USE_DEBUG_MALLOC 0
#define USE_ALIGNED_MALLOC 0
#define USE_FAST_MALLOC 0
#define USE_MY_DBGMALLOC 0

//
// Enable malloc debug
//
#if USE_DEBUG_MALLOC==1
//#include <malloc.h>
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif // USE_DEBUG_MALLOC


// [ malloc debug ] 

#define ALIGN_BYTES 8

#if USE_ALIGNED_MALLOC==1
#ifndef _CRTDBG_MAP_ALLOC
void * operator new(size_t size)
{
        void* p = _aligned_malloc(size, ALIGN_BYTES);
//        printf("Allocated: %x\n", p);
        return p;
};

void operator delete(void * mem)
{
//    printf("Delete %x\n", mem);
    _aligned_free(mem);
};

void * operator new[](size_t size)
{
        void* p = _aligned_malloc(size, ALIGN_BYTES);
//        printf("Allocated: %x\n", p);
        return p;
};
void operator delete[](void * mem)
{
//    printf("Delete %x\n", mem);
    _aligned_free(mem);
};
#endif // aligned malloc
#endif // USE_ALIGNED_MALLOC

// Use fast malloc
#if USE_FAST_MALLOC==1
void * operator new(size_t size)
{
        void* p = dlmalloc(size);
//        printf("Allocated: %x\n", p);
        return p;
};

void operator delete(void * mem)
{
//    printf("Delete %x\n", mem);
    dlfree(mem);
};

void * operator new[](size_t size)
{
        void* p = dlmalloc(size);
//        printf("Allocated: %x\n", p);
        return p;
};
void operator delete[](void * mem)
{
//    printf("Delete %x\n", mem);
    dlfree(mem);
};
#endif // USE_FAST_MALLOC==1


// My debug alloc
#if USE_MY_DBGMALLOC==1

struct blockinfo
{
    blockinfo() {}
    blockinfo(unsigned size_, int time_) : size(size_), time(time_) {}

    unsigned size;
    int time;
};

typedef std::map<void*, blockinfo> t_allocations;
static t_allocations allocations;
static bool collect_allocations=false;
static int current_time;

void start_leak_collection() 
{ 
    collect_allocations = true; 
    allocations.clear();
    current_time = 0;
}
void stop_leak_collection()
{
    collect_allocations = false;
}

void dump_leaks()
{
    t_allocations::iterator f = allocations.begin();
    t_allocations::iterator l = allocations.end();
    printf("Dumping leaks (%d) %s\n", current_time, f==l ? "No leaks" : "");
    for (; f !=l; ++f)
    {
        printf("0x%x\tsize: %d\ttime: %d\n", f->first, f->second.size, f->second.time);
    }
}

void * operator new(size_t size)
{
    void* p = malloc(size);
    if (collect_allocations)
    {
        collect_allocations = false;
        allocations[p] = blockinfo(size, ++current_time);
        collect_allocations = true;
    }
    return p;
};

void operator delete(void * mem)
{
    free(mem);
    if (collect_allocations)
    {
        collect_allocations = false;
        allocations.erase(mem);
        collect_allocations = true;
    }
};

void * operator new[](size_t size)
{
    void* p = malloc(size);
    if (collect_allocations)
    {
        collect_allocations = false;
        allocations[p] = blockinfo(size, ++current_time);
        collect_allocations = true;
    }
    return p;
};
void operator delete[](void * mem)
{
    free(mem);
    if (collect_allocations)
    {
        collect_allocations = false;
        allocations.erase(mem);
        collect_allocations = true;
    }
};
#else
void start_leak_collection() {}
void stop_leak_collection() {}
void dump_leaks() {}
#endif // USE_MY_DBGMALLOC==1

