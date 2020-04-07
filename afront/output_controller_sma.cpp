
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
#include "output_controller_sma.h"

#include <sm/smwriteopener.h>


OutputControllerSMA::OutputControllerSMA(const char *filename) {

	is_sma = (strstr(filename, ".sma") != NULL);

	SMwriteOpener opener;
	opener.set_file_name(filename);
	writer = opener.open();
	nverts=0;
	lasttri[0]=lasttri[1]=lasttri[2] = -1;
	finished[0]=finished[1]=finished[2] = false;
};


void OutputControllerSMA::AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary)
{
    if (writer) {
		float v[3] = {p[0],p[1],p[2]};
		writer->write_vertex(v);
    }
	nverts++;

    if (child)
	child->AddVertex(index, p, n, boundary);
}

void OutputControllerSMA::AddTriangle(int index, int v1, int v2, int v3)
{
    if (v1==v2 || v1==v3 || v2==v3)
	cerr<<"bogus triangle"<<endl;

	FlushLastTri();
	lasttri[0] = v1;
	lasttri[1] = v2;
	lasttri[2] = v3;


    if (child)
		child->AddTriangle(index, v1, v2, v3);
}


void OutputControllerSMA::FinalizeVertex(int index) {

	if (is_sma) {
 		FlushLastTri();
		writer->write_finalized(index);
	} else {
		for (int i=0; i<3; i++) {
			if (lasttri[i] == index)
				finished[i] = true;
		}
	}

    if (child) child->FinalizeVertex(index);
}

void OutputControllerSMA::Finish() {
	FlushLastTri();
	if (writer)
		writer->close();

	if (child)
		child->Finish();
}


void OutputControllerSMA::FlushLastTri() {

    if (writer && lasttri[0]>=0) {
		// unfinished: index+1, finished:index-vcount

		writer->write_triangle(lasttri, finished);

		// for (int i=0; i<3; i++) {

		// 	(*fout) << " ";
		// 	if (finished[i])
		// 		(*fout) << (lasttri[i]-nverts);
		// 	else
		// 		(*fout) << (lasttri[i]+1);
		// }
		// (*fout)<<endl;
    }


	lasttri[0] = lasttri[1] = lasttri[2] = -1;
	finished[0] = finished[1] = finished[2] = false;
}
