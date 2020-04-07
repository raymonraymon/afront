
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


#ifndef __MYNR_H
#define __MYNR_H

/*
 * Stuff from numerical recipes I use
 */

/*------------------ NR CODE ---------------*/
/*
 * Stuff from numerical recipes
 */
extern "C" float ran1f(long *idum);
extern "C" long ran1(long *idum);
extern "C" double gasdev(long* idum);



/*-------------- My extensions -------------*/
GTB_BEGIN_NAMESPACE

void nrseed(long v);
long nrran1();
float nrran1f();
double nrgasdev();

GTB_END_NAMESPACE

#endif // __MYNR_H
