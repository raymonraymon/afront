
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
#include <gtb/io/io.hpp>
#include <stdlib.h>
#endif // WIN32
using namespace gtb;


int main(int argc, char *argv[])
{
	if ((2 != argc) && (3 != argc)) {
		fprintf(stderr, "Usage: %s <FILE> [SUFFIX]\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	char fname[PATH_MAX];
	if (2 == argc) {
		get_file_base_name(argv[1], fname, PATH_MAX);
	} else if (3 == argc) {
		get_file_base_name(argv[1], argv[2], fname, PATH_MAX);
	} else {
		fprintf(stderr, "internal error\n");
		exit(EXIT_FAILURE);
	}
  	printf("%s\n", fname);
	return 0;
}
