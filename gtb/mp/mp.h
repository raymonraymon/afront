
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


#ifndef __GTB_META_PROGRAMMING_H
#define __GTB_META_PROGRAMMING_H

#include <gtb/common.hpp>
#include <gtb/mp/mp_common.h>

GTB_BEGIN_NAMESPACE
GTB_MP_BEGIN_NAMESPACE

#include "mp_factorial.inl"
#include "mp_power.inl"
#include "mp_choose.inl"
#include "mp_int2type.inl"
#include "mp_ctassert.h"

GTB_MP_END_NAMESPACE
GTB_END_NAMESPACE

#endif
