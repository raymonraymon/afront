
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


#ifndef __RUNINTHREAD_H
#define __RUNINTHREAD_H

#include "Thread.h"

BEGIN_THREADLIB_NAMESPACE

/*
 * This class is a helper class to run a method
 * in a thread.
 *
 * Template argument: the class that we run the thread in
 * ARG - the argument type to pass to the method
 *       Support for single argument only
 */

template <class T, class ARG>
class CRunInThread
{
public:
    typedef void (T::* mem_fun_type)(ARG);

    CRunInThread(T* p_class, mem_fun_type p_mem, ARG arg, int priority = THREAD_PRIORITY_NORMAL);
    virtual ~CRunInThread();

    void WaitForThread();

protected:
    T* _p_class;
    mem_fun_type _p_mem;
    ARG& _arg;
    HANDLE h_thread;
    static void ThreadEntry(CRunInThread* that);
    void runit();
};

/*
 * A similar class to the above, just that
 * the caller can continue doing its business 
 * without worrying about any of the stuff that is going on
 * including cleanup.
 *
 * The caller supplies a function to call and a pointer to delete
 * when the called function finishes execution
 */

template <class T, class ARG, class CLEANUP> 
class CRunInThreadAndForget : public CRunInThread<T, ARG>
{
public:
    CRunInThreadAndForget(T* p_class, mem_fun_type p_mem, ARG arg, CLEANUP* p_object);

protected:
    CLEANUP* _p_object;
    static void LocalThreadEntry(CRunInThreadAndForget* that);
    void runit();
};

template <class T, class ARG, class CLEANUP>
inline void RunInThreadAndForget(T* p_class, void (T::* p_mem)(ARG), ARG arg, CLEANUP* p_object)
{
    CRunInThreadAndForget<T, ARG, CLEANUP>* o_thread = 
        new CRunInThreadAndForget<T, ARG, CLEANUP>(p_class, p_mem, arg, p_object);
}


END_THREADLIB_NAMESPACE

#include "RunInThread.inl"

#endif // __RUNINTHREAD_H
