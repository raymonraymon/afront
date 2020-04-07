
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



#ifndef __OUTPUT_CONTROLLER_EDGEFLIPPER_H
#define __OUTPUT_CONTROLLER_EDGEFLIPPER_H

#include "triangulator.h"



class EdgeFlipperVertex {
    public:
    
    bool boundary;
    bool onfront;

    std::set<int> incidentTriangles;
    std::set<int> incidentVerticesLocal;
    std::set<int> incidentVerticesGlobal;
    int index;  // global index

    Point3 position;
    Vector3 normal;
};


class EdgeFlipperTriangle {
    public:
    int vi[3];  // local indices
    int index;  // global index
};





class OutputControllerEdgeFlipper : public OutputController
{
public:
    virtual ~OutputControllerEdgeFlipper() { }

    virtual void AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary);
    virtual void AddTriangle(int index, int v1, int v2, int v3);
    virtual void FinalizeVertex(int index);
    virtual void Finish();


    void Draw();

protected:

 
    vector<EdgeFlipperVertex> verts;
    vector<EdgeFlipperTriangle> tris;

    std::map<int,int> global_to_local_vi;

    vector<int> available_verts;
    vector<int> available_tris;
    int AllocVert();
    int AllocTri();
    void FreeVert(int v);
    void FreeTri(int t);

    void FinishVertex(int v);
    void FinishTriangle(int t);

    bool VertexRemovable(int vi) const;
    bool TriangleFlippable(int ti) const;

    bool TryFlip(int v1, int v2);
    bool ShouldFlip(int va, int vb, int vc, int vd) const;

    void DoEdgeFlips(const vector<int> &rv);
};


#endif
