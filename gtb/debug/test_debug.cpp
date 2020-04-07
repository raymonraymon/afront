
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
#include <stdlib.h>
#endif // WIN32
using namespace gtb;


void count(int n)
{
	debug d("counting");
	d.print("i: ");
	for (int i = 0; i < n; i++) {
		d.update("i: %d", i);
	}
	d.print("\n");
}


void foo()
{
	debug d("foo");
	count(100000);
}


int main(int, char *[])
{
	debug d("main");
	foo();
	return EXIT_SUCCESS;
}
