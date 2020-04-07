
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


#ifndef __HIGHRESTIMER_H
#define __HIGHRESTIMER_H

#include <gtb/common.hpp>
#include <stdio.h>
#include <time.h>
#include <string>

GTB_BEGIN_NAMESPACE

/*
 * Name: timer.h
 * Author: Shachar Fleishman
 * Description:
 *   High-resolution timer class
 *
 * Change
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 *
 */


#ifdef WIN32
/*
 * Timer class to help measure timing of program execution
 * time is measured in units of the machine.
 */
class CHighResTimer
{
public:

    typedef __int64 t_clock_ticks;

    CHighResTimer(bool auto_start = false) : m_elapsed(0) { if (auto_start) start(); }
    void start() { m_start = get_ticks(); }
    void stop() { m_elapsed += get_ticks() - m_start; }
    t_clock_ticks pause() {t_clock_ticks cur = get_ticks(); m_elapsed += cur - m_start; m_start = cur; return cur;}
    t_clock_ticks sample() {t_clock_ticks cur = get_ticks(); m_elapsed = cur - m_start; return m_elapsed;}
    void reset() { m_elapsed = 0; }
    operator t_clock_ticks() const { return m_elapsed; }
    int secs() const;
    int milisecs() const;
    void print() const { printf("%d.%d", secs(), milisecs()); }

    static t_clock_ticks get_ticks()
    {
        t_clock_ticks now;
        QueryPerformanceCounter((LARGE_INTEGER*)&now);
        return now;
	}

    static t_clock_ticks ticks_per_second()
    {
        t_clock_ticks tps;
        QueryPerformanceFrequency((LARGE_INTEGER*)&tps);
        return tps;
    }

protected:
    t_clock_ticks m_start;			// Start tick
    t_clock_ticks m_elapsed;		// Time elappsed (1/1000 sec)

};

inline int CHighResTimer::secs() const
{
    return m_elapsed / ticks_per_second(); 
}

inline int CHighResTimer::milisecs() const 
{ 
    t_clock_ticks tps = ticks_per_second();
    return (m_elapsed % tps) * 1000 / tps;
}

#else // NON WIN32

class CHighResTimer
{
public:

    typedef int t_clock_ticks;

    CHighResTimer(bool auto_start = false) : m_elapsed(0) { if (auto_start) start(); }
    void start() { m_start = get_ticks(); }
    void stop() { m_elapsed += get_ticks() - m_start; }
    t_clock_ticks pause() {t_clock_ticks cur = get_ticks(); m_elapsed += cur - m_start; m_start = cur; return cur;}
    t_clock_ticks sample() {t_clock_ticks cur = get_ticks(); m_elapsed = cur - m_start; return m_elapsed;}
    void reset() { m_elapsed = 0; }
    operator t_clock_ticks() const { return m_elapsed; }
    int secs() const { return 0; }
    int milisecs() const { return 0; }
    void print() const { printf("%d.%d", secs(), milisecs()); }

    static t_clock_ticks get_ticks()
    {
        return 0;
	}
protected:
    t_clock_ticks m_start;			// Start tick
    t_clock_ticks m_elapsed;		// Time elappsed (1/1000 sec)

};

#endif // WIN32

class HighResAutoStop
{
public:
    HighResAutoStop(CHighResTimer& timer) : _timer(timer) { _timer.start(); }
    ~HighResAutoStop() { _timer.stop(); }
protected:
    CHighResTimer& _timer;
};

/*
 * A helper class for accumulating time over several iterations
 * Upon construction of this class, the referenced timer begins counting
 * Upon destruction of this class, the reference timer stops
 * So, it makes timing exception safe.
 */
class CHRTimerAccumulator
{
public:
    CHRTimerAccumulator(CHighResTimer& t) : _t(t) { _t.start(); }
    ~CHRTimerAccumulator() { _t.stop(); }

private:
    CHighResTimer& _t;
};

GTB_END_NAMESPACE

#endif //__HIGHRESTIMER_H
