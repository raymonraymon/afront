
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


#ifndef GTB_TIMER_INCLUDED
#define GTB_TIMER_INCLUDED

#include <gtb/common.hpp>
//#include <stdlib.h>

#ifdef WIN32
//#include <winsock2.h>
#else
#include <sys/time.h>
#endif

GTB_BEGIN_NAMESPACE


class timer {
public:
	timer();
	void start();
	void stop();
	double elapsed_seconds() const;
	double elapsed_milliseconds() const;
protected:
	struct timeval m_t1;
	struct timeval m_t2;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/util/timer.ipp>
#endif

#endif // GTB_TIMER_INCLUDED
