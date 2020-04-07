
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
#include "output_controller_hhm.h"

void OutputControllerHHM::AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary)
{
    if (fout) {
		// 0 index not allowed? increase all indiced by 1 for .m files
		(*fout) << "Vertex  "<<(index+1)<<" "<<p[0]<<" "<<p[1]<<" "<<p[2];
		if (vkey)
			(*fout) <<" "<<vkey;

		(*fout)<<endl;
    }

    if (child)
	child->AddVertex(index, p, n, boundary);
}

void OutputControllerHHM::AddTriangle(int index, int v1, int v2, int v3)
{

    if (v1==v2 || v1==v3 || v2==v3)
		cerr<<"bogus triangle: "<<v1<<" "<<v2<<" "<<v3<<endl;
    if (fout) {
	// 0 index not allowed? increase all indiced by 1 for .m files
	(*fout) << "Face  "<<(index+1)<<" "<<(v1+1)<<" "<<(v2+1)<<" "<<(v3+1)<<endl;
    }

    if (child)
	child->AddTriangle(index, v1, v2, v3);
}

