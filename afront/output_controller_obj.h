
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


#ifndef __OUTPUT_CONTROLLER_OBJ_H
#define __OUTPUT_CONTROLLER_OBJ_H

#include "triangulator.h"
#include <fstream>

class OutputControllerOBJ : public OutputController
{
public:
    OutputControllerOBJ(const char *filename):
      fout(filename?new std::ofstream(filename):NULL) {};
    virtual ~OutputControllerOBJ() { delete fout; };
    virtual void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary);
    virtual void AddTriangle(int index, int v1, int v2, int v3);
    virtual void FinalizeVertex(int index) { if (child) child->FinalizeVertex(index); }
    virtual void Finish() { if (child) child->Finish(); }


protected:
    std::ofstream *fout;
};

#endif
