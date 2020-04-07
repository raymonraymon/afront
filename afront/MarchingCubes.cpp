
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



#include "common.h"
#include "parallel.h"
#include "guidance.h"
#include "triangulator.h"
#include "triangulate_iso.h"


#include <lib/mlslib/NR/nr.h>	// for zbrent?!?!?



class LineEval {

    public:
    LineEval(const Volume &v, const Point3 &_p0, const Point3 &_p1, real_type iso) :
		volume(v), p0(_p0), p1(_p1), iso_value(iso) { last_cell[0] = last_cell[1] = last_cell[2] = -1; }

    real_type operator()(real_type t) {

		Point3 x = p0 + t*(p1-p0);

		int this_cell[3];
		double xl[3];	// local coordinate

		volume.LocalInfo(x, this_cell, xl);
	    
		if (this_cell[0]!=last_cell[0] || this_cell[1]!=last_cell[1] || this_cell[2]!=last_cell[2]) {
	      
			volume.Gather(this_cell, nbrs);
			last_cell[0]=this_cell[0];
			last_cell[1]=this_cell[1];
			last_cell[2]=this_cell[2];
		}

		double e = volume.spline.Eval(volume.GetAspect(0),
									  volume.GetAspect(1),
									  volume.GetAspect(2),
									  nbrs, xl);

		return (real_type)(e - iso_value);
    }

        
    const Volume &volume;
    int last_cell[3]; 
    double nbrs[4][4][4];

    real_type iso_value;
    Point3 p0;
    Point3 p1;
};






#define	MARCH_NOP	0
#define MARCH_INVERT	1
#define MARCH_ROTATE_X	2
#define MARCH_ROTATE_Y	3
#define MARCH_ROTATE_Z	4



int mc_transforms[48][7] =
{

	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},

	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},

	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},

	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},

	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},

	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		},


	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},

	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},

	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},

	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_Y,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},

	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},

	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_NOP,		MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_NOP,		MARCH_INVERT,		},
	{ MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_Z,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_ROTATE_X,	MARCH_INVERT,		},
};


int mc_transform_index(int index, int trans)
{
	int ret = 0;

	switch (trans)
		{
		case MARCH_NOP:
			ret = index;
			break;

		case MARCH_INVERT:
			ret = ~index;
			break;

		case MARCH_ROTATE_X:
			if (index & (1<<0)) ret |= (1<<3);
			if (index & (1<<1)) ret |= (1<<2);
			if (index & (1<<2)) ret |= (1<<6);
			if (index & (1<<3)) ret |= (1<<7);
			if (index & (1<<4)) ret |= (1<<0);
			if (index & (1<<5)) ret |= (1<<1);
			if (index & (1<<6)) ret |= (1<<5);
			if (index & (1<<7)) ret |= (1<<4);
			break;

		case MARCH_ROTATE_Y:
			if (index & (1<<0)) ret |= (1<<4);
			if (index & (1<<1)) ret |= (1<<0);
			if (index & (1<<2)) ret |= (1<<3);
			if (index & (1<<3)) ret |= (1<<7);
			if (index & (1<<4)) ret |= (1<<5);
			if (index & (1<<5)) ret |= (1<<1);
			if (index & (1<<6)) ret |= (1<<2);
			if (index & (1<<7)) ret |= (1<<6);
			break;

		case MARCH_ROTATE_Z:
			if (index & (1<<0)) ret |= (1<<1);
			if (index & (1<<1)) ret |= (1<<2);
			if (index & (1<<2)) ret |= (1<<3);
			if (index & (1<<3)) ret |= (1<<0);
			if (index & (1<<4)) ret |= (1<<5);
			if (index & (1<<5)) ret |= (1<<6);
			if (index & (1<<6)) ret |= (1<<7);
			if (index & (1<<7)) ret |= (1<<4);
			break;

		default:
			cerr<<"mc_transform_index: unknown transform "<<trans<<endl;
			ret = index;
		}

	return ret;
}


void mc_transform_tri(int *tris, int num, int trans)
{
	int i, tmp;
	int xrot[12] = { 8, 5, 0, 4, 11, 9, 1, 3, 10, 6, 2, 7 };
	int yrot[12] = { 5, 9, 6, 1, 0, 8, 10, 2, 4, 11, 7, 3 };
	int zrot[12] = { 3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10 };


	switch (trans)
		{
		case MARCH_NOP:
			break;

		case MARCH_INVERT:
			for (i=0; i<(num*3); i+=3)
				{
					tmp = tris[i];
					tris[i] = tris[i+1];
					tris[i+1] = tmp;
				}
			break;

		case MARCH_ROTATE_X:
			for (i=0; i<(num*3); i++)
				tris[i] = xrot[tris[i]];
			break;

		case MARCH_ROTATE_Y:
			for (i=0; i<(num*3); i++)
				tris[i] = yrot[tris[i]];
			break;

		case MARCH_ROTATE_Z:
			for (i=0; i<(num*3); i++)
				tris[i] = zrot[tris[i]];
			break;

		default:
			cerr<<"transform_tri: unknown transform "<<trans<<endl;
		}
}


bool mc_testmatch(int index, int trans, bool extra, int *tris, int *ntris)
{
	int i;

	if ((index == 0) || (index == 255))
		{
			(*ntris) = 0;
			return true;
		}

	for (i=0; i<7; i++)
		index = mc_transform_index(index, mc_transforms[trans][i]);

	(*ntris) = 0;


	if (true)
		{
			if (index == 0)
				{
					(*ntris) = 0;
					return true;
				}

			else if (index == (1<<4))
				{
					(*ntris) = 1;

					tris[0] = 4;
					tris[1] = 11;
					tris[2] = 8;
					return true;
				}

			else if (index == ((1<<4) | (1<<5)))
				{
					(*ntris) = 2;

					tris[0] = 4;
					tris[1] = 11;
					tris[2] = 9;

					tris[3] = 4;
					tris[4] = 9;
					tris[5] = 5;
					return true;
				}

			else if (index == ((1<<4) | (1<<6)))
				{
					(*ntris) = 2;

					tris[0] = 4;
					tris[1] = 11;
					tris[2] = 8;

					tris[3] = 6;
					tris[4] = 9;
					tris[5] = 10;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<5)))
				{
					(*ntris) = 3;

					tris[0] = 4;
					tris[1] = 8;
					tris[2] = 9;

					tris[3] = 4;
					tris[4] = 9;
					tris[5] = 3;

					tris[6] = 3;
					tris[7] = 9;
					tris[8] = 1;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<4) | (1<<5)))
				{
					(*ntris) = 2;

					tris[0] = 11;
					tris[1] = 1;
					tris[2] = 3;

					tris[3] = 11;
					tris[4] = 9;
					tris[5] = 1;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<5) | (1<<7)))
				{
					(*ntris) = 4;

					tris[0] = 4;
					tris[1] = 8;
					tris[2] = 9;

					tris[3] = 4;
					tris[4] = 9;
					tris[5] = 3;

					tris[6] = 3;
					tris[7] = 9;
					tris[8] = 1;

					tris[9 ] = 11;
					tris[10] = 7;
					tris[11] = 10;
					return true;
				}

			else if (index == ((1<<1) | (1<<3) | (1<<4) | (1<<6)))
				{
					(*ntris) = 4;

					tris[0] = 11;
					tris[1] = 8;
					tris[2] = 4;

					tris[3] = 0;
					tris[4] = 5;
					tris[5] = 1;

					tris[6] = 2;
					tris[7] = 7;
					tris[8] = 3;

					tris[9 ] = 6;
					tris[10] = 9;
					tris[11] = 10;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<3) | (1<<4)))
				{
					(*ntris) = 4;

					tris[0] = 11;
					tris[1] = 2;
					tris[2] = 7;

					tris[3] = 11;
					tris[4] = 8;
					tris[5] = 2;

					tris[6] = 2;
					tris[7] = 8;
					tris[8] = 1;

					tris[9 ] = 1;
					tris[10] = 8;
					tris[11] = 5;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<3) | (1<<5)))
				{
					(*ntris) = 4;

					tris[0] = 7;
					tris[1] = 4;
					tris[2] = 8;

					tris[3] = 7;
					tris[4] = 8;
					tris[5] = 1;

					tris[6] = 7;
					tris[7] = 1;
					tris[8] = 2;

					tris[9] = 8;
					tris[10] = 9;
					tris[11] = 1;


					return true;
				}

			else if (index == ((1<<2) | (1<<4)))
				{
					(*ntris) = 2;

					tris[0] = 11;
					tris[1] = 8;
					tris[2] = 4;

					tris[3] = 6;
					tris[4] = 2;
					tris[5] = 1;
					return true;
				}

			else if (index == ((1<<2) | (1<<4) | (1<<5)))
				{
					(*ntris) = 3;

					tris[0] = 6;
					tris[1] = 2;
					tris[2] = 1;

					tris[3] = 11;
					tris[4] = 9;
					tris[5] = 4;

					tris[6] = 4;
					tris[7] = 9;
					tris[8] = 5;
					return true;
				}

			else if (index == ((1<<2) | (1<<5) | (1<<7)))
				{
					(*ntris) = 3;

					tris[0] = 6;
					tris[1] = 2;
					tris[2] = 1;

					tris[3] = 10;
					tris[4] = 11;
					tris[5] = 7;

					tris[6] = 8;
					tris[7] = 9;
					tris[8] = 5;
					return true;
				}

			else if (index == ((1<<1) | (1<<2) | (1<<4) | (1<<7)))
				{
					(*ntris) = 4;

					tris[0] = 2;
					tris[1] = 0;
					tris[2] = 6;

					tris[3] = 6;
					tris[4] = 0;
					tris[5] = 5;

					tris[6] = 7;
					tris[7] = 10;
					tris[8] = 4;

					tris[9 ] = 4;
					tris[10] = 10;
					tris[11] = 8;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<2) | (1<<4)))
				{
					(*ntris) = 4;

					tris[0] = 11;
					tris[1] = 8;
					tris[2] = 3;

					tris[3] = 3;
					tris[4] = 8;
					tris[5] = 6;

					tris[6] = 6;
					tris[7] = 8;
					tris[8] = 5;

					tris[9 ] = 3;
					tris[10] = 6;
					tris[11] = 2;
					return true;
				}

			//		else
			//			return false;
		}

	if (extra)
		{
			if (index == ((1<<0) | (1<<1) | (1<<3) | (1<<4) | (1<<5)))
				{
					(*ntris) = 3;

					tris[0] = 7;
					tris[1] = 11;
					tris[2] = 9;

					tris[3] = 7;
					tris[4] = 9;
					tris[5] = 1;

					tris[6] = 7;
					tris[7] = 1;
					tris[8] = 2;
					return true;
				}

			else if (index == ((1<<0) | (1<<2) | (1<<3) | (1<<4) | (1<<5)))
				{
					(*ntris) = 5;

					tris[0] = 1;
					tris[1] = 6;
					tris[2] = 0;

					tris[3] = 6;
					tris[4] = 11;
					tris[5] = 0;

					tris[6] = 0;
					tris[7] = 11;
					tris[8] = 9;

					tris[9 ] = 0;
					tris[10] = 9;
					tris[11] = 5;

					tris[12] = 7;
					tris[13] = 11;
					tris[14] = 6;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<2) | (1<<3) | (1<<4) | (1<<5)))
				{
					(*ntris) = 2;

					tris[0] = 7;
					tris[1] = 11;
					tris[2] = 9;

					tris[3] = 7;
					tris[4] = 9;
					tris[5] = 6;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<3) | (1<<4) | (1<<6)))
				{
					(*ntris) = 5;

					tris[0] = 11;
					tris[1] = 8;
					tris[2] = 7;

					tris[3] = 7;
					tris[4] = 8;
					tris[5] = 2;

					tris[6] = 2;
					tris[7] = 8;
					tris[8] = 5;

					tris[9 ] = 2;
					tris[10] = 5;
					tris[11] = 1;

					tris[12] = 9;
					tris[13] = 10;
					tris[14] = 6;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<3) | (1<<4) | (1<<5) | (1<<6)))
				{
					(*ntris) = 4;

					tris[0] = 7;
					tris[1] = 11;
					tris[2] = 2;

					tris[3] = 2;
					tris[4] = 11;
					tris[5] = 1;

					tris[6] = 1;
					tris[7] = 11;
					tris[8] = 10;

					tris[9 ] = 1;
					tris[10] = 10;
					tris[11] = 6;
					return true;
				}

			else if (index == ((1<<0) | (1<<2) | (1<<3) | (1<<4) | (1<<5) | (1<<6)))
				{
					(*ntris) = 2;

					tris[0] = 7;
					tris[1] = 11;
					tris[2] = 10;

					tris[3] = 0;
					tris[4] = 1;
					tris[5] = 5;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<2) | (1<<3) | (1<<4) | (1<<5) | (1<<6)))
				{
					(*ntris) = 1;

					tris[0] = 7;
					tris[1] = 11;
					tris[2] = 10;
					return true;
				}

			else if (index == ((1<<0) | (1<<1) | (1<<2) | (1<<3) | (1<<4) | (1<<5) | (1<<6) | (1<<7)))
				{
					(*ntris) = 0;
					return true;
				}

		}

	return false;
}


int64 mc_xyz2gindex(const Volume &v, int x, int y, int z) {
    return ((int64)z*(v.GetDim(0)+1)*(v.GetDim(1)+1) +
			(int64)y*(v.GetDim(0)+1) +
			(int64)x);
}

int mc_xyz2bindex(const int bsize[3], int x, int y, int z) {
    return (z*(bsize[0]+1)*(bsize[1]+1) +
			y*(bsize[0]+1) +
			x);
}

void mc_bindex2xyz(const int bsize[3], int i, int &x, int &y, int &z) {
	x = i % (bsize[0]+1);
	i /= (bsize[0]+1);
	
	y = i % (bsize[1]+1);
	i /= (bsize[1]+1);
	
	z = i;
}

int& mc_vindex(const Volume &v, HASH_MAP64 &gindices,
			   HASH_MAP32 &bindices, const int block_min[3], const int block_max[3], const int bsize[3],
			   int x, int y, int z, int c) {

	bool onbound = (x<=block_min[0] || x>=block_max[0]-1 || 
                  y<=block_min[1] || y>=block_max[1]-1 || 
                  z<=block_min[2] || z>=block_max[2]-1);

    // use global indices
    if (onbound) {
		return gindices[mc_xyz2gindex(v,x,y,z)*3 + c];
    } else {
		return bindices[mc_xyz2bindex(bsize, x-block_min[0], y-block_min[1], z-block_min[2])*3 + c];
    }
}




void MarchingCubesParallelVerts(int nt, int id, CSObject &cs, 
								const Volume &v, TriangulatorController &tc, 
								real_type isovalue, bool rootfinding, int mc_block, 
								int &num_verts,
								HASH_MAP64 &gindices, HASH_MAP32 &bindices,
								vector<int> *vertex_types) {

	vector<int> my_index;
	vector<Point3> my_points;
	vector<int> my_types;
	

	int block_min[3];
	int block_max[3];
	v.GetBlockBox(block_min, block_max);

	if (mc_block >= 0) {
		for (int i=0; i<3; i++) {
			block_max[i] = std::min(block_max[i]+1, v.GetDim(i));
		}
	}

	real_type bmin, bmax;
	v.GetBlockMinMax(bmin, bmax);


    int min_indices[3];
    int max_indices[3];

    if (mc_block >= 0) {
		v.SetActiveBlock(mc_block);
		v.GetBlockBox(min_indices, max_indices);
		
		for (int i=0; i<3; i++) {
			if (max_indices[i] == v.GetDim(i))
				max_indices[i]--;
		}
    } else {
		for (int i=0; i<3; i++) {
			min_indices[i] = 0;
			max_indices[i] = v.GetDim(i)-1;
		}
    }



	bool block_spans = (bmin<=isovalue && bmax>=isovalue);

	int bsize[3] = { block_max[0]-block_min[0],
					 block_max[1]-block_min[1],
					 block_max[2]-block_min[2] };

	bool hadgriderror = false;

	int i=0;
	for (ZCurve curve(bsize[0], bsize[1], bsize[2]); !curve.Done(); ++curve,++i) {

		int bx, by, bz;
		curve.Get(bx,by,bz);

		if (i%nt != id)
			continue;


		if (!block_spans) {
			bool bound = ((bx<2 || bx>bsize[0]-3) ||
						  (by<2 || by>bsize[1]-3) ||
						  (bz<2 || bz>bsize[2]-3));
			if (!bound)
				continue;
		}


		int vx = block_min[0] + bx;
		int vy = block_min[1] + by;
		int vz = block_min[2] + bz;

		real_type vcellvals[2][2][2];
		int vcell[3] = {vx,vy,vz};

		bool ccv = v.CellCornerValues(vcell, vcellvals, isovalue);
		{

			if (vcellvals[0][0][0] == isovalue ||
				vcellvals[0][0][1] == isovalue ||
				vcellvals[0][1][0] == isovalue ||
				vcellvals[0][1][1] == isovalue ||
				vcellvals[1][0][0] == isovalue ||
				vcellvals[1][0][1] == isovalue ||
				vcellvals[1][1][0] == isovalue ||
				vcellvals[1][1][1] == isovalue) {

				if (!hadgriderror) {
					cerr<<"grid point exactly equal to isovalue: "<<isovalue<<endl;
					cerr<<"you'll probably want to jitter the isovalue by a small amount (~1e-5)"<<endl;
					hadgriderror = true;
				}
			}

			// x-edge
			if ((vx < max_indices[0]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[1][0][0]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[1][0][0] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*(vx+0),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*(vx+1),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				my_index.push_back(mc_xyz2bindex(bsize, vx-block_min[0],vy-block_min[1],vz-block_min[2])*3+0);
				my_points.push_back(Point3(v.GetAspect(0)*(vx+s),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation());

				if (vertex_types) {
					int type = 0;
					if (vy==min_indices[1])    type |= MC_VERT_TYPE_YMIN;
					if (vy==max_indices[1])    type |= MC_VERT_TYPE_YMAX;
					if (vz==min_indices[2])    type |= MC_VERT_TYPE_ZMIN;
					if (vz==max_indices[2])    type |= MC_VERT_TYPE_ZMAX;
					my_types.push_back(type);
				}
			}

			// y-edge
			if ((vy < max_indices[1]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[0][1][0]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[0][1][0] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+0),v.GetAspect(2)*vz) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+1),v.GetAspect(2)*vz) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				my_index.push_back(mc_xyz2bindex(bsize, vx-block_min[0],vy-block_min[1],vz-block_min[2])*3+1);
				my_points.push_back(Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+s),v.GetAspect(2)*vz) + v.GetTranslation());

				if (vertex_types) {
					int type = 0;
					if (vx==min_indices[0])    type |= MC_VERT_TYPE_XMIN;
					if (vx==max_indices[0])    type |= MC_VERT_TYPE_XMAX;
					if (vz==min_indices[2])    type |= MC_VERT_TYPE_ZMIN;
					if (vz==max_indices[2])    type |= MC_VERT_TYPE_ZMAX;
					my_types.push_back(type);
				}
			}

			// z-edge
			if ((vz < max_indices[2]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[0][0][1]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[0][0][1] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+0)) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+1)) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				my_index.push_back(mc_xyz2bindex(bsize, vx-block_min[0],vy-block_min[1],vz-block_min[2])*3+2);
				my_points.push_back(Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+s)) + v.GetTranslation());

				if (vertex_types) {
					int type = 0;
					if (vx==min_indices[0])    type |= MC_VERT_TYPE_XMIN;
					if (vx==max_indices[0])    type |= MC_VERT_TYPE_XMAX;
					if (vy==min_indices[1])    type |= MC_VERT_TYPE_YMIN;
					if (vy==max_indices[1])    type |= MC_VERT_TYPE_YMAX;
					my_types.push_back(type);
				}
			}
		}
	}

	// copy everything back in one chunk
	cs.enter();
	for (unsigned i=0; i<my_points.size(); i++) {
		
		int vx, vy, vz, vi, tmp;
		vi = my_index[i] % 3;
		mc_bindex2xyz(bsize, my_index[i]/3, vx, vy, vz);
		
		tc.AddVertex(num_verts, my_points[i], Vector3(0,0,0), false);
		mc_vindex(v, gindices, bindices, block_min, block_max, bsize, vx+block_min[0], vy+block_min[1], vz+block_min[2], vi) = num_verts;
		num_verts++;

		if (vertex_types)
			vertex_types->push_back(my_types[i]);
	}		 	
	cs.leave();
}



void MarchingCubesParallelTris(int nt, int id, CSObject &cs, 
							   const Volume &v, TriangulatorController &tc, 
							   real_type isovalue, bool rootfinding, int mc_block, 
							   int &num_tris,
							   HASH_MAP64 &gindices, HASH_MAP32 &bindices,
							   vector<int64> *intersected_cells) {

	vector<int> my_tris;
	vector<int> my_cells;

	int block_min[3];
	int block_max[3];
	v.GetBlockBox(block_min, block_max);

	if (mc_block >= 0) {
		for (int i=0; i<3; i++) {
			block_max[i] = std::min(block_max[i]+1, v.GetDim(i));
		}
	}

	real_type bmin, bmax;
	v.GetBlockMinMax(bmin, bmax);


    int min_indices[3];
    int max_indices[3];

    if (mc_block >= 0) {
		v.SetActiveBlock(mc_block);
		v.GetBlockBox(min_indices, max_indices);
		
		for (int i=0; i<3; i++) {
			if (max_indices[i] == v.GetDim(i))
				max_indices[i]--;
		}
    } else {
		for (int i=0; i<3; i++) {
			min_indices[i] = 0;
			max_indices[i] = v.GetDim(i)-1;
		}
    }



	bool block_spans = (bmin<=isovalue && bmax>=isovalue);

	int bsize[3] = { block_max[0]-block_min[0],
					 block_max[1]-block_min[1],
					 block_max[2]-block_min[2] };


	int i=0;
	for (ZCurve curve(bsize[0], bsize[1], bsize[2]); !curve.Done(); ++curve,++i) {

		int bx, by, bz;
		curve.Get(bx,by,bz);

		if (i%nt != id)
			continue;


		if (!block_spans) {
			bool bound = ((bx<2 || bx>bsize[0]-3) ||
						  (by<2 || by>bsize[1]-3) ||
						  (bz<2 || bz>bsize[2]-3));
			if (!bound)
				continue;
		}


		int tx = block_min[0] + bx - 1;
		int ty = block_min[1] + by - 1;
		int tz = block_min[2] + bz - 1;

		if (tx < min_indices[0]) continue;
		if (ty < min_indices[1]) continue;
		if (tz < min_indices[2]) continue;


		real_type tcellvals[2][2][2];
		int tcell[3] = {tx,ty,tz};

		v.CellCornerValues(tcell, tcellvals, isovalue);
		{

			int index = 0;
			if (tcellvals[0][0][0]> isovalue)	index |= 1<<0;
			if (tcellvals[1][0][0]> isovalue)	index |= 1<<1;
			if (tcellvals[1][1][0]> isovalue)	index |= 1<<2;
			if (tcellvals[0][1][0]> isovalue)	index |= 1<<3;
			if (tcellvals[0][0][1]> isovalue)	index |= 1<<4;
			if (tcellvals[1][0][1]> isovalue)	index |= 1<<5;
			if (tcellvals[1][1][1]> isovalue)	index |= 1<<6;
			if (tcellvals[0][1][1]> isovalue)	index |= 1<<7;


			// first check the extra cases with no inversion
			bool vol_march_found = false;

			int n_tris;
			int tris[5*3];
	

			int trans;
			for (trans=0; trans<24; trans++) {
				if (mc_testmatch(index, trans, true, tris, &n_tris)) {
					vol_march_found=true; break;
				}
			}
			    
			// check the reglar cases, and do the inversions
			if (!vol_march_found) {
				for (trans = 0; trans<48; trans++) {
					if (mc_testmatch(index, trans, false, tris, &n_tris)) {
						vol_march_found=true; break;
					}
				}
			}

			if (!vol_march_found)
				cerr<<"didn't find match for cube!"<<endl;

			if (intersected_cells!=NULL && n_tris>0) {
				int64 cell = tz;
				cell *= v.GetDim(1);
				cell += ty;
				cell *= v.GetDim(0);
				cell += tx;

				my_cells.push_back(cell);
			}

			// rotate all the tris back
			for (int t=6; t>=0; t--) {
				mc_transform_tri(tris, n_tris, mc_transforms[trans][t]);
			}

			for (int t=0; t<(n_tris*3); t++) {
				switch (tris[t]) {
				case 0:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,0); break;
				case 1:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz  ,1); break;
				case 2:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz  ,0); break;
				case 3:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,1); break;
				case 4:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,2); break;
				case 5:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz  ,2); break;
				case 6:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty+1,tz  ,2); break;
				case 7:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz  ,2); break;
				case 8:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz+1,0); break;
				case 9:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz+1,1); break;
				case 10:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz+1,0); break;
				case 11:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz+1,1); break;
				}
			}

			for (int t=0; t<n_tris*3; t++) {
				my_tris.push_back(tris[t]);
			}

		}
	}

	cs.enter();
	for (unsigned t=0; t<my_tris.size(); t+=3) {
		tc.AddTriangle(num_tris++, my_tris[t+0], my_tris[t+1], my_tris[t+2]);
	}
	if (intersected_cells) {
		for (unsigned i=0; i<my_cells.size(); i++) 
			intersected_cells->push_back(my_cells[i]);
	}
	cs.leave();
}


// can't do streaming in parallel!
void MarchingCubesStreaming(const Volume &v, TriangulatorController &tc, 
							real_type isovalue, bool rootfinding, int mc_block, 
							HASH_MAP64 &gindices, HASH_MAP32 &bindices,
							int &num_verts, int &num_tris,
							vector<int64> *intersected_cells, vector<int> *vertex_types) {

	int block_min[3];
	int block_max[3];
	v.GetBlockBox(block_min, block_max);

	if (mc_block >= 0) {
		for (int i=0; i<3; i++) {
			block_max[i] = std::min(block_max[i]+1, v.GetDim(i));
		}
	}

	real_type bmin, bmax;
	v.GetBlockMinMax(bmin, bmax);


    int min_indices[3];
    int max_indices[3];

    if (mc_block >= 0) {
		v.SetActiveBlock(mc_block);
		v.GetBlockBox(min_indices, max_indices);
		
		for (int i=0; i<3; i++) {
			if (max_indices[i] == v.GetDim(i))
				max_indices[i]--;
		}
    } else {
		for (int i=0; i<3; i++) {
			min_indices[i] = 0;
			max_indices[i] = v.GetDim(i)-1;
		}
    }



	bool block_spans = (bmin<=isovalue && bmax>=isovalue);

	int bsize[3] = { block_max[0]-block_min[0],
					 block_max[1]-block_min[1],
					 block_max[2]-block_min[2] };

	bool hadgriderror = false;

	int i=0;
	for (ZCurve curve(bsize[0], bsize[1], bsize[2]); !curve.Done(); ++curve,++i) {

		int bx, by, bz;
		curve.Get(bx,by,bz);

		if (!block_spans) {
			bool bound = ((bx<2 || bx>bsize[0]-3) ||
						  (by<2 || by>bsize[1]-3) ||
						  (bz<2 || bz>bsize[2]-3));
			if (!bound)
				continue;
		}


		int vx = block_min[0] + bx;
		int vy = block_min[1] + by;
		int vz = block_min[2] + bz;

		real_type vcellvals[2][2][2];
		int vcell[3] = {vx,vy,vz};

		bool ccv = v.CellCornerValues(vcell, vcellvals, isovalue);
		{

			if (vcellvals[0][0][0] == isovalue ||
				vcellvals[0][0][1] == isovalue ||
				vcellvals[0][1][0] == isovalue ||
				vcellvals[0][1][1] == isovalue ||
				vcellvals[1][0][0] == isovalue ||
				vcellvals[1][0][1] == isovalue ||
				vcellvals[1][1][0] == isovalue ||
				vcellvals[1][1][1] == isovalue) {

				if (!hadgriderror) {
					cerr<<"grid point exactly equal to isovalue: "<<isovalue<<endl;
					cerr<<"you'll probably want to jitter the isovalue by a small amount (~1e-5)"<<endl;
					hadgriderror = true;
				}
			}

			// x-edge
			if ((vx < max_indices[0]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[1][0][0]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[1][0][0] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*(vx+0),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*(vx+1),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				tc.AddVertex(num_verts, Point3(v.GetAspect(0)*(vx+s),v.GetAspect(1)*vy,v.GetAspect(2)*vz) + v.GetTranslation(), Vector3(0,0,0), false);
				mc_vindex(v, gindices, bindices, block_min, block_max, bsize, vx, vy, vz, 0) = num_verts;
				num_verts++;


				if (vertex_types) {
					int type = 0;
					if (vy==min_indices[1])    type |= MC_VERT_TYPE_YMIN;
					if (vy==max_indices[1])    type |= MC_VERT_TYPE_YMAX;
					if (vz==min_indices[2])    type |= MC_VERT_TYPE_ZMIN;
					if (vz==max_indices[2])    type |= MC_VERT_TYPE_ZMAX;
					vertex_types->push_back(type);
				}
			}

			// y-edge
			if ((vy < max_indices[1]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[0][1][0]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[0][1][0] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+0),v.GetAspect(2)*vz) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+1),v.GetAspect(2)*vz) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				tc.AddVertex(num_verts, Point3(v.GetAspect(0)*vx,v.GetAspect(1)*(vy+s),v.GetAspect(2)*vz) + v.GetTranslation(), Vector3(0,0,0), false);
				mc_vindex(v, gindices, bindices, block_min, block_max, bsize, vx, vy, vz, 1) = num_verts;
				num_verts++;

				if (vertex_types) {
					int type = 0;
					if (vx==min_indices[0])    type |= MC_VERT_TYPE_XMIN;
					if (vx==max_indices[0])    type |= MC_VERT_TYPE_XMAX;
					if (vz==min_indices[2])    type |= MC_VERT_TYPE_ZMIN;
					if (vz==max_indices[2])    type |= MC_VERT_TYPE_ZMAX;
					vertex_types->push_back(type);
				}
			}

			// z-edge
			if ((vz < max_indices[2]) && (vcellvals[0][0][0]<isovalue) ^ (vcellvals[0][0][1]<isovalue)) {
				real_type s = (isovalue - vcellvals[0][0][0]) / (vcellvals[0][0][1] - vcellvals[0][0][0]);

				if (rootfinding) {
					Point3 p1 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+0)) + v.GetTranslation();
					Point3 p2 = Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+1)) + v.GetTranslation();
					LineEval le(v, p1, p2, isovalue);

					extern real_type step_tol;
					extern real_type value_tol;
					real_type step = step_tol;

					for (int i=0; i<5; i++) {
						s = NR::zbrent(le, (real_type)0, (real_type)1, step);
						step *= 0.05;
						if (fabs(le(s)) < value_tol)
							break;
					}
				}

				tc.AddVertex(num_verts, Point3(v.GetAspect(0)*vx,v.GetAspect(1)*vy,v.GetAspect(2)*(vz+s)) + v.GetTranslation(), Vector3(0,0,0), false);
				mc_vindex(v, gindices, bindices, block_min, block_max, bsize, vx, vy, vz, 2) = num_verts;
				num_verts++;

				if (vertex_types) {
					int type = 0;
					if (vx==min_indices[0])    type |= MC_VERT_TYPE_XMIN;
					if (vx==max_indices[0])    type |= MC_VERT_TYPE_XMAX;
					if (vy==min_indices[1])    type |= MC_VERT_TYPE_YMIN;
					if (vy==max_indices[1])    type |= MC_VERT_TYPE_YMAX;
					vertex_types->push_back(type);
				}
			}
		}



		int tx = vx-1;
		int ty = vy-1;
		int tz = vz-1;

		if (tx < min_indices[0]) continue;
		if (ty < min_indices[1]) continue;
		if (tz < min_indices[2]) continue;


		real_type tcellvals[2][2][2];
		int tcell[3] = {tx,ty,tz};

		v.CellCornerValues(tcell, tcellvals, isovalue);
		{

			int index = 0;
			if (tcellvals[0][0][0]> isovalue)	index |= 1<<0;
			if (tcellvals[1][0][0]> isovalue)	index |= 1<<1;
			if (tcellvals[1][1][0]> isovalue)	index |= 1<<2;
			if (tcellvals[0][1][0]> isovalue)	index |= 1<<3;
			if (tcellvals[0][0][1]> isovalue)	index |= 1<<4;
			if (tcellvals[1][0][1]> isovalue)	index |= 1<<5;
			if (tcellvals[1][1][1]> isovalue)	index |= 1<<6;
			if (tcellvals[0][1][1]> isovalue)	index |= 1<<7;


			// first check the extra cases with no inversion
			bool vol_march_found = false;

			int n_tris;
			int tris[5*3];
	

			int trans;
			for (trans=0; trans<24; trans++) {
				if (mc_testmatch(index, trans, true, tris, &n_tris)) {
					vol_march_found=true; break;
				}
			}
			    
			// check the reglar cases, and do the inversions
			if (!vol_march_found) {
				for (trans = 0; trans<48; trans++) {
					if (mc_testmatch(index, trans, false, tris, &n_tris)) {
						vol_march_found=true; break;
					}
				}
			}

			if (!vol_march_found)
				cerr<<"didn't find match for cube!"<<endl;

			if (intersected_cells!=NULL && n_tris>0) {
				int64 cell = tz;
				cell *= v.GetDim(1);
				cell += ty;
				cell *= v.GetDim(0);
				cell += tx;

				intersected_cells->push_back(cell);
			}

			// rotate all the tris back
			for (int t=6; t>=0; t--) {
				mc_transform_tri(tris, n_tris, mc_transforms[trans][t]);
			}

			bool final_x = false;
			bool final_y = false;
			bool final_z = false;

			for (int t=0; t<(n_tris*3); t++) {
				switch (tris[t]) {
				case 0:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,0); final_x=true; break;
				case 1:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz  ,1); break;
				case 2:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz  ,0); break;
				case 3:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,1); final_y=true; break;
				case 4:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,2); final_z=true; break;
				case 5:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz  ,2); break;
				case 6:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty+1,tz  ,2); break;
				case 7:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz  ,2); break;
				case 8:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz+1,0); break;
				case 9:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx+1,ty  ,tz+1,1); break;
				case 10:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty+1,tz+1,0); break;
				case 11:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz+1,1); break;
				}
			}

			for (int t=0; t<n_tris*3; t+=3) {
				tc.AddTriangle(num_tris++, tris[t+0], tris[t+1], tris[t+2]);
			}

			if (final_x)
				tc.FinalizeVertex(mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,0));
			if (final_y)
				tc.FinalizeVertex(mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,1));
			if (final_z)
				tc.FinalizeVertex(mc_vindex(v,gindices,bindices,block_min,block_max,bsize, tx  ,ty  ,tz  ,2));

		}
	}
}




void MarchingCubes(const Volume &v, TriangulatorController &tc, real_type isovalue, bool rootfinding, int mc_block, vector<int64> *intersected_cells, vector<int> *vertex_types) {


	extern bool stream_mc;

    int num_verts=0;
    int num_tris =0;
	int &num_verts_r = num_verts;
	int &num_tris_r = num_tris;

    HASH_MAP64 gindices;
    HASH_MAP64 &gindices_r = gindices;


    // find all the vertices (on the edges) 
    for (int block = ((mc_block>=0) ? mc_block : 0);
		 block < ((mc_block>=0) ? mc_block+1 : v.NumBlocks());
		 block++) {

		HASH_MAP32 bindices;
		HASH_MAP32 &bindices_r = bindices;
	
		int block_min[3];
		int block_max[3];
		v.SetActiveBlock(block);
		v.GetBlockBox(block_min, block_max);

		if (mc_block >= 0) {
			for (int i=0; i<3; i++) {
				block_max[i] = std::min(block_max[i]+1, v.GetDim(i));
			}
		}

		real_type bmin, bmax;
		v.GetBlockMinMax(bmin, bmax);


		bool block_spans = (bmin<=isovalue && bmax>=isovalue);

		int bsize[3] = { block_max[0]-block_min[0],
						 block_max[1]-block_min[1],
						 block_max[2]-block_min[2] };


		if (stream_mc) {
			MarchingCubesStreaming(v, tc, isovalue, rootfinding, mc_block,
								   gindices_r, bindices_r,
								   num_verts_r, num_tris_r,
								   intersected_cells, vertex_types);
		} else {
			ParallelExecutor(idealNumThreads, MarchingCubesParallelVerts, 
							 v, tc, isovalue, rootfinding, mc_block, num_verts_r, gindices_r, bindices_r, vertex_types);
			ParallelExecutor(idealNumThreads, MarchingCubesParallelTris, 
							 v, tc, isovalue, rootfinding, mc_block, num_tris_r, gindices_r, bindices_r, intersected_cells);
		}
    }    
}




void MarchingCubesStats(const Volume &v, TriangulatorController *tc, real_type isovalue,
						real_type &cellcount, real_type &wcellcount,
						real_type &tricount, real_type &wtricount,
						real_type &area, real_type &warea) {


	cellcount = 0;
	wcellcount = 0;
	tricount = 0;
	wtricount = 0;
	area = 0;
	warea = 0;


    int min_indices[3];
    int max_indices[3];


	bool rootfinding=false;
	int mc_block=-1;

	int num_verts = 0;
	int num_tris = 0;
 
    if (mc_block >= 0) {
		v.SetActiveBlock(mc_block);
		v.GetBlockBox(min_indices, max_indices);

		for (int i=0; i<3; i++) {
			if (max_indices[i] == v.GetDim(i))
				max_indices[i]--;
		}
    } else {
		for (int i=0; i<3; i++) {
			min_indices[i] = 0;
			max_indices[i] = v.GetDim(i)-1;
		}
    }


    // find all the vertices (on the edges) 
    for (int block = ((mc_block>=0) ? mc_block : 0);
		 block < ((mc_block>=0) ? mc_block+1 : v.NumBlocks());
		 block++) {

		int block_min[3];
		int block_max[3];
		v.SetActiveBlock(block);
		v.GetBlockBox(block_min, block_max);

		if (mc_block >= 0) {
			for (int i=0; i<3; i++) {
				block_max[i] = std::min(block_max[i]+1, v.GetDim(i));
			}
		}

		real_type bmin, bmax;
		v.GetBlockMinMax(bmin, bmax);

		// this is not 100% safe for bspline interpolation!
		bool block_spans = (bmin<isovalue && bmax>isovalue);

		int bsize[3] = { block_max[0]-block_min[0],
						 block_max[1]-block_min[1],
						 block_max[2]-block_min[2] };

		//	bindices.resize((bsize[0]+1)*(bsize[1]+1)*(bsize[2]+1)*3);
		bool hadgriderror = false;

		for (ZCurve curve(bsize[0], bsize[1], bsize[2]); !curve.Done(); ++curve) {

			int bx, by, bz;
			curve.Get(bx,by,bz);

			if (!block_spans) {
				bool bound = ((bx<2 || bx>bsize[0]-3) ||
							  (by<2 || by>bsize[1]-3) ||
							  (bz<2 || bz>bsize[2]-3));
				if (!bound)
					continue;
			}


			int vx = block_min[0] + bx - 1;
			int vy = block_min[1] + by - 1;
			int vz = block_min[2] + bz - 1;

			real_type vcellvals[2][2][2];
			int vcell[3] = {vx,vy,vz};
			v.CellCornerValues(vcell, vcellvals, isovalue);
			{

				if (vcellvals[0][0][0] == isovalue ||
					vcellvals[0][0][1] == isovalue ||
					vcellvals[0][1][0] == isovalue ||
					vcellvals[0][1][1] == isovalue ||
					vcellvals[1][0][0] == isovalue ||
					vcellvals[1][0][1] == isovalue ||
					vcellvals[1][1][0] == isovalue ||
					vcellvals[1][1][1] == isovalue) {

					if (!hadgriderror) {
						cerr<<"grid point exactly equal to isovalue"<<endl;
						cerr<<"you'll probably want to jitter the isovalue by a small amount (~1e-5)"<<endl;
						hadgriderror = true;
					}
				}


				// check cell intersections
				real_type min_corner = vcellvals[0][0][0];
				real_type max_corner = vcellvals[0][0][0];
				for (int x=0; x<2; x++) {
					for (int y=0; y<2; y++) {
						for (int z=0; z<2; z++) {
							min_corner = std::min(min_corner, vcellvals[x][y][z]);
							max_corner = std::max(max_corner, vcellvals[x][y][z]);
						}
					}
				}			

				if (min_corner>=isovalue || max_corner<=isovalue)
					continue;


				Vector3 grad(0,0,0);
				if (min_corner<=isovalue && max_corner>=isovalue) {
					for (int y=0; y<2; y++) {
						for (int z=0; z<2; z++) {
							grad[0] += 0.25 * (vcellvals[1][y][z] - vcellvals[0][y][z]);
						}
					}

					for (int x=0; x<2; x++) {
						for (int z=0; z<2; z++) {
							grad[1] += 0.25 * (vcellvals[x][1][z] - vcellvals[x][0][z]);
						}
					}
			
					for (int x=0; x<2; x++) {
						for (int y=0; y<2; y++) {
							grad[2] += 0.25 * (vcellvals[x][y][1] - vcellvals[x][y][0]);
						}
					}
				}


				// reset these on every step!
				HASH_MAP64 gindices;
				HASH_MAP32 bindices;

				std::map<int,Point3> pmap;
				std::map<int,Vector3> gmap;
		


				// x-edges
				for (int y=0; y<2; y++) {
					for (int z=0; z<2; z++) {
						if ((vx < max_indices[0]) && (vcellvals[0][y][z]<isovalue) ^ (vcellvals[1][y][z]<isovalue)) {
							real_type s = (isovalue - vcellvals[0][y][z]) / (vcellvals[1][y][z] - vcellvals[0][y][z]);

							Point3 p(Point3(v.GetAspect(0)*(vx+s),
											v.GetAspect(1)*(vy+y),
											v.GetAspect(2)*(vz+z)) + v.GetTranslation());
							Vector3 n;
							v.GradientAtPoint(p, n);
							if (tc)
								tc->AddVertex(num_verts, p, n, false);
							mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx,vy+y,vz+z,0) = num_verts;
							pmap[num_verts] = p;
							gmap[num_verts] = n;
							num_verts++;
						}
					}
				}

				// y-edge
				for (int x=0; x<2; x++) {
					for (int z=0; z<2; z++) {
						if ((vy < max_indices[1]) && (vcellvals[x][0][z]<isovalue) ^ (vcellvals[x][1][z]<isovalue)) {
							real_type s = (isovalue - vcellvals[x][0][z]) / (vcellvals[x][1][z] - vcellvals[x][0][z]);

							Point3 p(Point3(v.GetAspect(0)*(vx+x),
											v.GetAspect(1)*(vy+s),
											v.GetAspect(2)*(vz+z)) + v.GetTranslation());
							Vector3 n;
							v.GradientAtPoint(p, n);
							if (tc)
								tc->AddVertex(num_verts, p, n, false);
							mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+x,vy,vz+z,1) = num_verts;
							pmap[num_verts] = p;
							gmap[num_verts] = n;
							num_verts++;
						}
					}
				}

				// z-edge
				for (int x=0; x<2; x++) {
					for (int y=0; y<2; y++) {
						if ((vz < max_indices[2]) && (vcellvals[x][y][0]<isovalue) ^ (vcellvals[x][y][1]<isovalue)) {
							real_type s = (isovalue - vcellvals[x][y][0]) / (vcellvals[x][y][1] - vcellvals[x][y][0]);

							Point3 p(Point3(v.GetAspect(0)*(vx+x),
											v.GetAspect(1)*(vy+y),
											v.GetAspect(2)*(vz+s)) + v.GetTranslation());
							Vector3 n;
							v.GradientAtPoint(p, n);
							if (tc)
								tc->AddVertex(num_verts, p, n, false);
							mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+x,vy+y,vz,2) = num_verts;
							pmap[num_verts] = p;
							gmap[num_verts] = n;
							num_verts++;

						}
					}
				}



				int index = 0;
				if (vcellvals[0][0][0]> isovalue)	index |= 1<<0;
				if (vcellvals[1][0][0]> isovalue)	index |= 1<<1;
				if (vcellvals[1][1][0]> isovalue)	index |= 1<<2;
				if (vcellvals[0][1][0]> isovalue)	index |= 1<<3;
				if (vcellvals[0][0][1]> isovalue)	index |= 1<<4;
				if (vcellvals[1][0][1]> isovalue)	index |= 1<<5;
				if (vcellvals[1][1][1]> isovalue)	index |= 1<<6;
				if (vcellvals[0][1][1]> isovalue)	index |= 1<<7;


				// first check the extra cases with no inversion
				bool vol_march_found = false;

				int n_tris;
				int tris[5*3];
	

				int trans;
				for (trans=0; trans<24; trans++) {
					if (mc_testmatch(index, trans, true, tris, &n_tris)) {
						vol_march_found=true; break;
					}
				}
			    
				// check the reglar cases, and do the inversions
				if (!vol_march_found) {
					for (trans = 0; trans<48; trans++) {
						if (mc_testmatch(index, trans, false, tris, &n_tris)) {
							vol_march_found=true; break;
						}
					}
				}

				if (!vol_march_found)
					cerr<<"didn't find match for cube!"<<endl;


				// rotate all the tris back
				for (int t=6; t>=0; t--) {
					mc_transform_tri(tris, n_tris, mc_transforms[trans][t]);
				}

				for (int t=0; t<(n_tris*3); t++) {
					switch (tris[t]) {
					case 0:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy  ,vz  ,0); break;
					case 1:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+1,vy  ,vz  ,1); break;
					case 2:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy+1,vz  ,0); break;
					case 3:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy  ,vz  ,1); break;
					case 4:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy  ,vz  ,2); break;
					case 5:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+1,vy  ,vz  ,2); break;
					case 6:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+1,vy+1,vz  ,2); break;
					case 7:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy+1,vz  ,2); break;
					case 8:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy  ,vz+1,0); break;
					case 9:	tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx+1,vy  ,vz+1,1); break;
					case 10:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy+1,vz+1,0); break;
					case 11:tris[t] = mc_vindex(v,gindices,bindices,block_min,block_max,bsize, vx  ,vy  ,vz+1,1); break;
					}
				}


				//		real_type gradmag = grad.length();
				real_type gradmag = max_corner-min_corner;

				real_type invgradmag = 1 / gradmag;
				if (gradmag == 0)
					invgradmag = 0;

				for (int t=0; t<n_tris*3; t+=3) {
					if (tc)
						tc->AddTriangle(num_tris++, tris[t+0], tris[t+1], tris[t+2]);

					Point3 p0 = pmap[tris[t+0]];
					Point3 p1 = pmap[tris[t+1]];
					Point3 p2 = pmap[tris[t+2]];

					real_type tarea = (p1-p0).cross(p2-p0).length() * 0.5;

					area += tarea;
					warea += tarea * invgradmag;

					tricount += 1;
					wtricount += 1 * invgradmag;
				}

				cellcount += 1;
				wcellcount += 1 * invgradmag;



			}
		}

    }    
}




