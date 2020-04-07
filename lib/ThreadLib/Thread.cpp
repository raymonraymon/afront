
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

#include "Thread.h"

#include <iostream>
#include <cassert>
using namespace std;

BEGIN_THREADLIB_NAMESPACE


#ifdef WIN32

typedef struct {
	void *(*thread_main)(void*);
	void *arg;
} winthreadstart_t;

DWORD WINAPI WinThreadStart(void *arg) {
	winthreadstart_t *wts = (winthreadstart_t*)arg;
	(*wts->thread_main)(wts->arg);
	delete wts;
	return 0;
}

Thread::Thread(void *(*where_to)(void *), void *parameter, int priority)
{
	 winthreadstart_t *wts = new winthreadstart_t;
	 wts->thread_main = where_to;
	 wts->arg = parameter;

	 _thread = CreateThread(NULL, 0, WinThreadStart, wts, 0, NULL);
//    SetThreadPriority(_thread, priority);
}

Thread::~Thread()
{
    CloseHandle(_thread);
}

void Thread::exit(void *status)
{
	ExitThread(*(DWORD*)status);
}

void Thread::join(void **)
{
    WaitForSingleObject(_thread, INFINITE);
}

void Thread::yield()
{
	SwitchToThread();
}

void Thread::sleep()
{
    Sleep(0);
}
#endif

#if defined(__gnu_linux__)
Thread::Thread(void *(*where_to)(void*), void *parameter, int)
{
    if (pthread_create(&_thread, NULL, where_to, parameter)) {
	cerr << "Thread::Thread(): pthread_create failed" << endl;
	assert(0);
    }
}

Thread::~Thread()
{
}

void Thread::exit(void *status)
{
    pthread_exit(status);
}

void Thread::join(void **status_p)
{
    if (pthread_join(_thread, status_p)) {
	cerr << "Thread::join(): pthread_join failed" << endl;
	assert(0);
    }
}

void Thread::yield()
{
    pthread_yield();
}

size_t Thread::self()
{
	return (size_t) pthread_self();
}
#endif

#if defined(__APPLE__)
Thread::Thread(void *(*where_to)(void*), void *parameter, int)
{
    if (pthread_create(&_thread, NULL, where_to, parameter)) {
	cerr << "Thread::Thread(): pthread_create failed" << endl;
	assert(0);
    }
}

Thread::~Thread()
{
}

void Thread::exit(void *status)
{
    pthread_exit(status);
}

void Thread::join(void **status_p)
{
    if (pthread_join(_thread, status_p)) {
	cerr << "Thread::join(): pthread_join failed" << endl;
	assert(0);
    }
}

void Thread::yield()
{
    sched_yield();
}
#endif

END_THREADLIB_NAMESPACE
