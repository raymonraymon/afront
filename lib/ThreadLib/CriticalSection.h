
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


#ifndef __CRITICALSECTION_H
#define __CRITICALSECTION_H

#include "threadsall.h"

/*
 * Critical section class
 * Sample usage:
 *
 * CSObject list_guard;
 * ...
 * in each thread do:
 *  {
 *     CS cs(list_guard); // Enter the critical section
 *     ...
 *  } // Releases the critical section
 */

BEGIN_THREADLIB_NAMESPACE

//
// A critical section object, used to guard
// an object. Must declare one for the object that
// is being guarded
//
class CSObject
{
public:
    CSObject();
    ~CSObject(); // NOTE: not virtual

    void enter();
    void leave();
private:
    CriticalSection _cs;
};


// spin lock
#ifndef __gnu_linux__

// win32 just use CSObjects
typedef CSObject SLObject;

#else

// otherwise use actual spins
class SLObject
{
public:
    SLObject() {
		pthread_spin_init(&lock,0);
	}
    ~SLObject() {
		pthread_spin_destroy(&lock);
	}

    void enter() {
		pthread_spin_lock(&lock);
	}
    void leave() {
		pthread_spin_unlock(&lock);
	}
private:
	pthread_spinlock_t lock;
};

#endif


//
// Auto release critical section class
//
class CS
{
public:
    CS(CSObject& csobject);
    ~CS();

    CSObject& _csobject;
};




END_THREADLIB_NAMESPACE


#endif // __CRITICALSECTION_H
