
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
#include "output_controller_obj.h"

void OutputControllerOBJ::AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary)
{
    if (fout) {
	// 0 index not allowed? increase all indiced by 1 for .m files
	(*fout) << "v "<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
	(*fout) << "vn "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl;
    }

    if (child)
	child->AddVertex(index, p, n, boundary);
}

void OutputControllerOBJ::AddTriangle(int index, int v1, int v2, int v3)
{
    if (fout) {
	// 0 index not allowed? increase all indiced by 1 for .m files
	(*fout) << "f "<<v1<<" "<<v2<<" "<<v3<<endl;
    }

    if (child)
	child->AddTriangle(index, v1, v2, v3);
}

