
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
#include "kdtree.h"

GTB_BEGIN_NAMESPACE

kdtree_exception_no_more_values _exception_no_more_values;

CHighResTimer _timer_push;
CHighResTimer _timer_pop;
CHighResTimer _timer_insert;
CHighResTimer _timer_makeheap;
CHighResTimer _timer_ordered_init_pq;
CHighResTimer _timer_ordered_init_pq2;
CHighResTimer _timer_ordered_getnext;

GTB_END_NAMESPACE
