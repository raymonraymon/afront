
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



#include "stdafx.h"

#include "CriticalSection.h"


#include <cassert>
#include <iostream>
using namespace std;

BEGIN_THREADLIB_NAMESPACE

#ifdef WIN32
CSObject::CSObject()
{
    InitializeCriticalSection(&_cs);
}

CSObject::~CSObject()
{
    DeleteCriticalSection(&_cs);
}

void CSObject::enter()
{
    EnterCriticalSection(&_cs);
}

void CSObject::leave()
{
    LeaveCriticalSection(&_cs);
}
#endif

#if defined(__gnu_linux__) || defined(__APPLE__)
#include <errno.h>
CSObject::CSObject()
{
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE_NP);
    pthread_mutex_init(&_cs, &attr);
}

CSObject::~CSObject()
{
    pthread_mutex_destroy(&_cs);
}

void CSObject::enter()
{
	/*
    if (pthread_mutex_lock(&_cs)) 
    {
	cerr << "CSObject::enter() failed" << endl;
	assert(0);
    };
	*/

	while (1) {
		int res = pthread_mutex_trylock(&_cs);
		if (res == EBUSY) {
			continue; // try again
		} else if (res) {
			cerr << "CSObject::enter() failed" << endl;
			assert(0);
		} else {
			return;
		}
	}
	
}

void CSObject::leave()
{
    int error = 0;
    if (error=pthread_mutex_unlock(&_cs)) {
		cerr << "WHHHAAAT????  CSObject::leave() failed" << endl;
		if (error == EINVAL) {
			cerr << "invalid mutex" << endl;
		}
		if (error == EPERM) {
			cerr << "not owner" << endl;
		}
		cerr << "error: " << error << flush;
		assert(0);
	}
}
#endif

CS::CS(CSObject& csobject) :
    _csobject(csobject)
{
    _csobject.enter();
}

CS::~CS()
{
    _csobject.leave();
}

END_THREADLIB_NAMESPACE
