
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
 * file:	FLF_structs.h
 * description:	Saves some structures used in the output and creation
 *              of the feature loop files.
 *
 * author:	Joel Daniels and Linh Ha
 * date:	May 2006
 *
 */

#ifndef _FLF_STRUCTS_H_
#define _FLF_STRUCTS_H_

typedef struct {
  float x_,y_,z_;
  bool isCorner_;
} InfoPoint;

typedef struct {
  unsigned int index_;
  float x_,y_,z_;
} InfoIndex;

#endif
