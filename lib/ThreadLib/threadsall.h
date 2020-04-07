
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


#ifndef __THREADSALL_H
#define __THREADSALL_H

#if defined(__gnu_linux__) || defined(__APPLE__)
#include <pthread.h>
#include <semaphore.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <winbase.h>
extern "C"
WINBASEAPI
BOOL
WINAPI
SwitchToThread(
    VOID
    );

#endif

#define BEGIN_THREADLIB_NAMESPACE namespace thlib {
#define END_THREADLIB_NAMESPACE }

BEGIN_THREADLIB_NAMESPACE

#if defined(__gnu_linux__) || defined(__APPLE__)
typedef pthread_mutex_t CriticalSection;
typedef ::sem_t PrimitiveSemaphore;
typedef pthread_t PrimitiveThread;
#endif

#if defined(__APPLE__)
#define PTHREAD_MUTEX_RECURSIVE_NP PTHREAD_MUTEX_RECURSIVE
#endif

#ifdef WIN32
typedef CRITICAL_SECTION CriticalSection;
typedef HANDLE PrimitiveSemaphore;
typedef HANDLE PrimitiveThread;
#endif

END_THREADLIB_NAMESPACE

#endif // __THREADSALL_H
