
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
#include <gtb/error/error.hpp>
#include <stdio.h>
#include <stdlib.h>
#endif // WIN32
using namespace gtb;


int main(int argc, char *argv[])
{
	set_program_name(argv[0]);

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <[1-4]>\n", program_name());
		exit(EXIT_FAILURE);
	}
	int testno = atoi(argv[1]);

	switch (testno) {
	case 1:
	{
		const char *fname = "/non/existent/file";
		FILE *fp = fopen(fname, "rb");
		if (NULL == fp) {
			GTB_ERROR(fname);
		}
		printf("should not see this\n");
		break;
	}
	case 2:
		GTB_WARNING("division by zero\n");
		printf("should see this\n");
		break;
	case 3:
		REQUIRE(2 + 2 == 5);
		break;
	case 4:
		ENSURE(2 + 2 == 3);
		break;
	default:
		GTB_ERROR("invalid test number");
		break;
	}

	return EXIT_SUCCESS;
}
