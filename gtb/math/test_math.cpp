
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
#include <gtb/math/math.hpp>
#endif // WIN32

using namespace std;
using namespace gtb;


int main(int, char *[])
{
	printf("%f\n", M_PI / 4.0);
	printf("%f\n", deg_to_rad(45.0));
	printf("%f\n", rad_to_deg(M_PI / 4.0));
	printf("%d\n", abs(1));
	printf("%f\n", abs(-1.2));
	printf("%f\n", abs(0.0));
	printf("%d\n", max(1, 3));
	printf("%f\n", min(3.0, 1.0));
	return EXIT_SUCCESS;
}
