
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


#ifndef GTB_ERROR_INCLUDED
#define GTB_ERROR_INCLUDED

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

const char *program_name();
void set_program_name(const char *argv0);

void error(const char *file, int line, const char *fmt, ...);
void warning(const char *file, int line, const char *fmt, ...);

#define GTB_ERROR(args) error(__FILE__, __LINE__, args)
#define GTB_WARNING(args) warning(__FILE__, __LINE__, args)

#define REQUIRE(condition) \
do { \
	if (!(condition)) { \
		error(__FILE__, __LINE__, "precondition violation: %s\n", \
		      # condition); \
	} \
} while (0)

#define ENSURE(condition) \
do { \
	if (!(condition)) { \
		error(__FILE__, __LINE__, "postcondition violation: %s\n", \
		      # condition); \
	} \
} while (0)

GTB_END_NAMESPACE

#include <gtb/error/err.h>

#endif // GTB_ERROR_INCLUDED
