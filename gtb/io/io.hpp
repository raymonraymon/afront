
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


#ifndef GTB_IO_INCLUDED
#define GTB_IO_INCLUDED

#include <gtb/common.hpp>
#include <gtb/error/error.hpp>
#include <stdio.h>
#include <limits.h>

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

GTB_BEGIN_NAMESPACE


void read_unsigned(unsigned *i, FILE *fp);
void read_int(int *i, FILE *fp);
void read_float(float *f, FILE *fp);
void read_bool(bool *b, FILE *fp);

void write_unsigned(unsigned i, FILE *fp);
void write_int(int i, FILE *fp);
void write_float(float f, FILE *fp);
void write_bool(bool b, FILE *fp);

void write_double(double d, FILE *fp);
void read_double(double *d, FILE *fp);

void get_file_base_name(const char *file_name, char *base_name, unsigned size);
void get_file_base_name(const char *file_name,
			const char *suffix,
			char *base_name,
			unsigned size);
void get_file_extension(const char *file_name,
			char *extension,
			unsigned size);
FILE *xfopen(const char *file_name, const char *mode);
bool file_exists(const char *file_name);
bool dir_exists(const char *dir_name);

FILE* xfopenv(const char* file_name, int max_backup=5);

GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/io/io.ipp>
#endif

#endif // GTB_IO_INCLUDED
