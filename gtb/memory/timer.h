
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


#ifndef __TIMER_H
#define __TIMER_H

#include <gtb/common.hpp>
#include <stdio.h>
#include <time.h>
#include <string>

GTB_BEGIN_NAMESPACE

/*
 * Name: timer.h
 * Author: Shachar Fleishman
 * Description:
 *   simple timer class
 *
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */



/*
 * Timer class to help measure timing of program execution
 * time is measured in 1/1000 of a second.
 */
class CTimer
{
public:
#ifdef WIN32
	typedef DWORD t_clock_ticks;
#else
    typedef clock_t t_clock_ticks;
#endif

    CTimer(bool auto_start = false) : m_elapsed(0) { if (auto_start) start(); }
    void start() { m_start = get_ticks(); }
    void stop() { m_elapsed += get_ticks() - m_start; }
    int  pause() {int cur = get_ticks(); m_elapsed += cur - m_start; m_start = cur; return cur;}
    int  sample() {int cur = get_ticks(); m_elapsed = cur - m_start; return m_elapsed;}
    void reset() { m_elapsed = 0; }
    operator int() const { return m_elapsed; }
    int secs() const { return m_elapsed / 1000; }
    int milisecs() const { return m_elapsed % 1000; }
    void print() const { printf("%d.%d", secs(), milisecs()); }

    static t_clock_ticks get_ticks()
    {
#ifdef WIN32
	    // clock() is too high resolution and overlapps to fast
	    // although it may be more accurate..
	    return GetTickCount();
#else
	    return 1000 * clock() / CLOCKS_PER_SEC;
#endif
	}
protected:
    t_clock_ticks m_start;			// Start tick
    t_clock_ticks m_elapsed;		// Time elappsed (1/1000 sec)

};

/*
 * A helper timer class for debug and user info
 * does the following:
 *  print a user friendly message on screen at init
 *  print ellapsed time when deleted (optionally with the user message again)
 */
class CTimerOnce
{
public:
    CTimerOnce(const char* message, bool Reprint = false) :
        t(true),
        msg(message),
        repreint_message(Reprint)
    {
        printf("%s ", msg.c_str());
    }

    ~CTimerOnce()
    {
        t.stop();
        if (repreint_message)
        {
            printf("%s ", msg.c_str());
        }
        t.print();
        printf("\n");
    }
private:
    CTimer t;
    std::string msg;
    bool repreint_message;
};

/*
 * A helper class for accumulating time over several iterations
 * Upon construction of this class, the referenced timer begins counting
 * Upon destruction of this class, the reference timer stops
 * So, it makes timing exception safe.
 */
class CTimerAccumulator
{
public:
    CTimerAccumulator(CTimer& t) : _t(t) { _t.start(); }
    ~CTimerAccumulator() { _t.stop(); }

private:
    CTimer& _t;
};

GTB_END_NAMESPACE

#endif //__TIMER_H
