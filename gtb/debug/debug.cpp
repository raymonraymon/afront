
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
#ifndef WIN32
#include <gtb/debug/debug.hpp>
#include <stdarg.h>
#include <stdio.h>
#endif // WIN32


GTB_BEGIN_NAMESPACE

int debug::_indent = 0;


debug::debug(const char *s)
	: _s(s)
{
	fprintf(stderr, "\n");
	indent();
	fprintf(stderr, "%s {\n", _s);
	_indent++;
}


debug::~debug()
{
	_indent--;
	indent();
	fprintf(stderr, "} // %s\n\n", _s);
}


void debug::print(const char *fmt, ...) const
{
	va_list args;
	va_start(args, fmt);
	indent();
	vfprintf(stderr, fmt, args);
	va_end(args);
}


void debug::update(const char *fmt, ...) const
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "\r");
	indent();
	vfprintf(stderr, fmt, args);
	va_end(args);
}


void debug::indent() const
{
	for (int i = 0; i < _indent; i++) {
		fprintf(stderr, "    ");
	}
}


GTB_END_NAMESPACE
