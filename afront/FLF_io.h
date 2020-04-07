
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
 * file:	FLF_io.h
 * description:	Feature Loop File converter provides the file
 *		write and read functions that will convert our
 *		extracted features into a format that is necessary
 *		for John and Carlos' triangulation method.
 *
 * author:	Joel Daniels and Linh Ha
 * date:	May 2006
 *
 */

#ifndef _FLF_IO_H_
#define _FLF_IO_H_

// -- INCLUDES -- //
#include "FLF_structs.h"
#include <stdio.h>
#include <fstream>
#include <vector>

/**
 * These functions write the features extracted from the point cloud
 * as feature loops or read feature loops into the standard vectors.
 * The feature loops are ordered lists of indices with associated
 * normals at each point.  These indices will access the point list
 * also given to the file.
 */
void FLF_write( const char *filename, const std::vector< InfoPoint > &ptList, const std::vector< std::vector< InfoIndex > > &featList );
void FLF_read( const char *filename, std::vector< InfoPoint > &ptList, std::vector< std::vector< InfoIndex > > &featList );

#endif
