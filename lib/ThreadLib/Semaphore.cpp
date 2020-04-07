
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


#include "Sem.h"

#include <iostream>
#include <cassert>
using namespace std;

#ifdef WIN32
#endif

#ifdef __gnu_linux_
Semaphore::Semaphore(int initial_value)
{
    if (sem_init(&_semaphore, 0, initial_value)) {
	cerr << "Semaphore::Semaphore(): sem_init failed" << endl;
	assert(0);
    };
}

Semaphore::~Semaphore()
{
    if (sem_destroy(&_semaphore)) {
	cerr << "Semaphore::~Semaphore(): sem_destroy failed" << endl;
	assert(0);
    };
}

void Semaphore::v()
{
    if (sem_post(&_semaphore)) {
	cerr << "Semaphore::v(): sem_post failed" << endl;
	assert(0);
    }
}

void Semaphore::p()
{
    if (sem_wait(&_semaphore)) {
	cerr << "Semaphore::p(): sem_wait failed" << endl;
	assert(0);
    }
}
#endif

#ifdef WIN32
Semaphore::Semaphore(int initial_value)
{
    _semaphore = CreateSemaphore(NULL, initial_value, 20000000000, NULL);
}

Semaphore::~Semaphore()
{
    CloseHandle(_semaphore);
}

void Semaphore::v()
{
    ReleaseSemaphore(_semaphore);
}

void Semaphore::p()
{
    WaitForSingleObject(_semaphore, INFINITE);
}
#endif
