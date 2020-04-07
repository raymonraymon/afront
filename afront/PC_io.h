
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


/**
 *
 * file:	PC_io.h
 * description:	Point Cloud io operations.  read and write files
 *		that store point cloud information.
 *
 * author:	Joel Daniels and Linh Ha
 * date:	May 2006
 *
 */

#ifndef _PC_IO_H_
#define _PC_IO_H_

// -- INCLUDES -- //
#include "FLF_structs.h"
#include <stdio.h>
#include <fstream>
#include <vector>

/**
 * These functions write the points in a cloud.  It is just a point
 * list.
 */
void PC_write( const char *filename, const std::vector< InfoPoint > &ptList );
void PC_read( const char *filename, std::vector< InfoPoint > &ptList );

#endif
