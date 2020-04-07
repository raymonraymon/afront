
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
#include "edgeflipper.h"



void OutputControllerEdgeFlipper::AddVertex(int index, const Point3 &p, const Vector3 &n, bool boundary) {

    int vil = AllocVert();
    EdgeFlipperVertex &v = verts[vil];
    v.index = index;
    v.boundary = boundary;
    v.onfront = true;
    v.incidentTriangles.clear();
    v.incidentVerticesLocal.clear();
    v.incidentVerticesGlobal.clear();
    v.position = p;
    v.normal = n;

    global_to_local_vi[index] = vil;


    // vertex gets sent downstream right away
    child->AddVertex(index, p, n, boundary);
}


void OutputControllerEdgeFlipper::FinalizeVertex(int index) {
    int lindex = global_to_local_vi[index];
    verts[lindex].onfront = false;


    // figure out which verts are ready to remove
    // this stuff could probably be fancier
    std::set<int> prv;
    if (VertexRemovable(lindex))
		prv.insert(lindex);

    for (std::set<int>::iterator vi=verts[lindex].incidentVerticesLocal.begin();
		 vi!=verts[lindex].incidentVerticesLocal.end();
		 ++vi) {
	
		if (VertexRemovable(*vi)) {
			prv.insert(*vi);
		}
    }
    
    
    vector<int> remove_verts;
    for (std::set<int>::iterator i=prv.begin(); i!=prv.end(); ++i) {
		if (VertexRemovable(*i)) {
			remove_verts.push_back(*i);
		}
    }

	DoEdgeFlips(remove_verts);


    for (unsigned i=0; i<remove_verts.size(); i++) {
		if (!verts[remove_verts[i]].boundary)
			FinishVertex(remove_verts[i]);
    }
}

void OutputControllerEdgeFlipper::AddTriangle(int index, int v1, int v2, int v3) {

    int til = AllocTri();
    EdgeFlipperTriangle &t = tris[til];

    t.vi[0] = global_to_local_vi[v1];
    t.vi[1] = global_to_local_vi[v2];
    t.vi[2] = global_to_local_vi[v3];
    t.index = index;


    // add incident edges and the new triangle
    int gvi[3] = { v1, v2, v3 };
    for (int i=0; i<3; i++) {
		EdgeFlipperVertex &v = verts[t.vi[i]];
		v.incidentTriangles.insert(til);
		v.incidentVerticesLocal.insert(t.vi[(i+1)%3]);
		v.incidentVerticesLocal.insert(t.vi[(i+2)%3]);
		v.incidentVerticesGlobal.insert(gvi[(i+1)%3]);
		v.incidentVerticesGlobal.insert(gvi[(i+2)%3]);
    }

    vector<int> flipverts;
	flipverts.push_back(t.vi[0]);
	flipverts.push_back(t.vi[1]);
	flipverts.push_back(t.vi[2]);

	DoEdgeFlips(flipverts);
}


void OutputControllerEdgeFlipper::Finish() {

    for (unsigned t=0; t<tris.size(); t++) {
		bool used=true;
		for (unsigned i=0; i<available_tris.size(); i++) {
			if (available_tris[i] == t) {
				used = false;
				break;
			}
		}

		if (used) {
			
			if (verts[tris[t].vi[0]].index != verts[tris[t].vi[1]].index &&
				verts[tris[t].vi[1]].index != verts[tris[t].vi[2]].index &&
				verts[tris[t].vi[2]].index != verts[tris[t].vi[0]].index)
				
				child->AddTriangle(tris[t].index,
								   verts[tris[t].vi[0]].index,
								   verts[tris[t].vi[1]].index,
								   verts[tris[t].vi[2]].index);
		}
    }

	
    tris.clear();
    verts.clear();
    available_tris.clear();
    available_verts.clear();


    child->Finish();
}



real_type angle(Point3 &a, Point3 &b, Point3 &c) {
    return acossafe((a-b).normalized().dot((c-b).normalized()));
}
bool OutputControllerEdgeFlipper::ShouldFlip(int va, int vb, int vc, int vd) const {

    Point3 pts[4] = { verts[va].position, 
					  verts[vb].position, 
					  verts[vc].position, 
					  verts[vd].position };

    real_type curmin = std::min(std::min(angle(pts[0], pts[1], pts[2]),
										 std::min(angle(pts[1], pts[2], pts[0]),
												  angle(pts[2], pts[0], pts[1]))),
								std::min(angle(pts[3], pts[1], pts[2]),
										 std::min(angle(pts[1], pts[2], pts[3]),
												  angle(pts[2], pts[3], pts[1]))));

    real_type flipmin= std::min(std::min(angle(pts[0], pts[1], pts[3]),
										 std::min(angle(pts[1], pts[3], pts[0]),
												  angle(pts[3], pts[0], pts[1]))),
								std::min(angle(pts[3], pts[2], pts[0]),
										 std::min(angle(pts[2], pts[0], pts[3]),
												  angle(pts[0], pts[3], pts[2]))));


    if (flipmin < curmin) {
		return false;
    }


    Vector3 on1 = (pts[1]-pts[0]).cross(pts[2]-pts[0]);
    Vector3 on2 = (pts[2]-pts[3]).cross(pts[1]-pts[3]);

    Vector3 nn1 = (pts[3]-pts[1]).cross(pts[0]-pts[1]);
    Vector3 nn2 = (pts[0]-pts[2]).cross(pts[3]-pts[2]);

    real_type old_area = on1.length()+on2.length();
    real_type new_area = nn1.length()+nn2.length();


	if (new_area > 1.1*old_area)
	 	return false;

    // on1.normalize();
    // on2.normalize();
    // nn1.normalize();
    // nn2.normalize();

    // const real_type min_dot=0;
    // if (on1.dot(nn1)<min_dot || on1.dot(nn2)<min_dot)
	// return false;

    // if (on2.dot(nn1)<min_dot || on2.dot(nn2)<min_dot)
	// return false;


    return true;
}

bool OutputControllerEdgeFlipper::TryFlip(int v1, int v2) {


    int t1=-1;
    int t2=-1;

    int o1=-1;
    int o2=-1;


    if (verts[v1].incidentVerticesGlobal.size() <= 3) {
	return false;
    }

    if (verts[v2].incidentVerticesGlobal.size() <= 3) {
	return false;
    }


    for (std::set<int>::iterator ti=verts[v1].incidentTriangles.begin();
	 ti!=verts[v1].incidentTriangles.end();
	 ++ti) {

	if (*ti < 0)
	    continue;

	for (int i=0; i<3; i++) {
	    if (tris[*ti].vi[i] == v1) {
		if (tris[*ti].vi[(i+2)%3] == v2) {
		    t1 = *ti;
		    o1 = tris[*ti].vi[(i+1)%3];
		}
		if (tris[*ti].vi[(i+1)%3] == v2) {
		    t2 = *ti;
		    o2 = tris[*ti].vi[(i+2)%3];
		}
	    }
	}
    }


    if (t1<0 || t2<0) {
	return false;
    }


    int nbound=0;
    if (verts[v1].boundary) nbound++;
    if (verts[v2].boundary) nbound++;
    if (verts[o1].boundary) nbound++;
    if (verts[o2].boundary) nbound++;

	// make sure not to flip it if it's a boundary edge
    // if (nbound >= 1) {
	// return false;
    // }
	if (verts[v1].boundary && verts[v2].boundary)
		return false;


    if (!ShouldFlip(o1, v2, v1, o2))
	return false;


    // remove the old stuff
    std::set<int>::iterator i;
    int vs[2] = { v1, v2 };
    for (int v=0; v<2; v++) {

	i = verts[vs[v]].incidentVerticesLocal.find(vs[1-v]);
	if (i == verts[vs[v]].incidentVerticesLocal.end()) {
	    cerr<<"flip: local vert not found!"<<endl;
	} else {
	    verts[vs[v]].incidentVerticesLocal.erase(i);
	}

	i = verts[vs[v]].incidentVerticesGlobal.find(verts[vs[1-v]].index);
	if (i == verts[vs[v]].incidentVerticesGlobal.end()) {
	    cerr<<"flip: global vert not found!"<<endl;
	} else {
	    verts[vs[v]].incidentVerticesGlobal.erase(i);
	}
    }

    i = verts[v1].incidentTriangles.find(t2);
    if (i == verts[v1].incidentTriangles.end()) {
	cerr<<"flip: v1 tri not found!"<<endl;
    } else {
	verts[v1].incidentTriangles.erase(i);
    }

    i = verts[v2].incidentTriangles.find(t1);
    if (i == verts[v2].incidentTriangles.end()) {
	cerr<<"flip: v2 tri not found!"<<endl;
    } else {
	verts[v2].incidentTriangles.erase(i);
    }


    // change where things are pointing
    bool found;
    found=false;
    for (int i=0; i<3; i++) {
	if (tris[t1].vi[i] == v2) {
	    tris[t1].vi[i] = o2;
	    found=true;
	    break;
	}
    }
    if (!found)
	cerr<<"wtf"<<endl;

    found=false;
    for (int i=0; i<3; i++) {
	if (tris[t2].vi[i] == v1) {
	    tris[t2].vi[i] = o1;
	    found=true;
	    break;
	}
    }
    if (!found)
	cerr<<"wtf"<<endl;


    verts[o1].incidentVerticesLocal.insert(o2);
    verts[o2].incidentVerticesLocal.insert(o1);

    verts[o1].incidentVerticesGlobal.insert(verts[o2].index);
    verts[o2].incidentVerticesGlobal.insert(verts[o1].index);

    verts[o1].incidentTriangles.insert(t2);
    verts[o2].incidentTriangles.insert(t1);


    return true;
}


void OutputControllerEdgeFlipper::DoEdgeFlips(const vector<int> &rv) {

    std::set< std::pair<int,int> > toflip;

    while (1) {

	for (unsigned v=0; v<rv.size(); v++) {

	    for (std::set<int>::iterator ti=verts[rv[v]].incidentTriangles.begin();
		 ti!=verts[rv[v]].incidentTriangles.end();
		 ++ti) {

		if (*ti > 0) {
		    for (int i=0; i<3; i++) {
			std::pair<int,int> edge(tris[*ti].vi[(i+0)%3], tris[*ti].vi[(i+1)%3]);
			if (edge.second < edge.first)
			    std::swap(edge.first, edge.second);
			toflip.insert(edge);
		    }
		}

	    }
	}


	bool didflip = false;
	while (!toflip.empty()) {
	    std::pair<int,int> edge = *toflip.begin();
	    toflip.erase(toflip.begin());
	    
	    didflip |= TryFlip(edge.first, edge.second);
	}

	if (!didflip) break;
    }


}



void OutputControllerEdgeFlipper::Draw() {


    for (unsigned t=0; t<tris.size(); t++) {
	bool inuse=true;

	for (unsigned i=0; i<available_tris.size(); i++) {
	    if (t == available_tris[i]) {
		inuse = false;
		break;
	    }
	}
	if (!inuse) continue;

	int vof = 0;
	for (int i=0; i<3; i++) {
	    if (verts[tris[t].vi[i]].onfront) {
		vof++;
	    }
	}

	vof = 0;

	Point3 color;
	switch (vof) {
	case 0:    color = Point3(0.5,0.5,0.5);    break;
	case 1:    color = Point3(0,0,1);    break;
	case 2:    color = Point3(0,1,0);    break;
	case 3:    color = Point3(1,0,0);    break;
	}

	Point3 centroid = Point3(0,0,0) +
	    (real_type)0.333333 * (verts[tris[t].vi[0]].position - Point3(0,0,0)) +
	    (real_type)0.333333 * (verts[tris[t].vi[1]].position - Point3(0,0,0)) +
	    (real_type)0.333333 * (verts[tris[t].vi[2]].position - Point3(0,0,0));

	real_type shrink = (real_type)1.0;
	DbgPoly::add(centroid + shrink * (verts[tris[t].vi[0]].position-centroid), color,
		     centroid + shrink * (verts[tris[t].vi[1]].position-centroid), color,
		     centroid + shrink * (verts[tris[t].vi[2]].position-centroid), color);


	vector<Point3> line(4);
	line[0] = centroid + shrink * (verts[tris[t].vi[0]].position-centroid);
	line[1] = centroid + shrink * (verts[tris[t].vi[1]].position-centroid);
	line[2] = centroid + shrink * (verts[tris[t].vi[2]].position-centroid);
	line[3] = centroid + shrink * (verts[tris[t].vi[0]].position-centroid);
	DbgPLines::add(line);
	

	/*
	DbgPoly::add(verts[tris[t].vi[0]].position, color,
		     verts[tris[t].vi[1]].position, color,
		     verts[tris[t].vi[2]].position, color);
	*/

    }
#if 0
    for (unsigned v=0; v<verts.size(); v++) {
	bool inuse=true;
	for (unsigned i=0; i<available_verts.size(); i++) {
	    if (v == available_verts[i]) {
		inuse = false;
		break;
	    }
	}
	if (!inuse) continue;


	if (verts[v].onfront && verts[v].boundary)
	    DbgPoints::add(verts[v].position, 1, 0, 0);

	if (verts[v].onfront && !verts[v].boundary)
	    DbgPoints::add(verts[v].position, 0.5, 0, 0);

	if (!verts[v].onfront && verts[v].boundary)
	    DbgPoints::add(verts[v].position, 0, 1, 0);

	if (!verts[v].onfront && !verts[v].boundary)
	    DbgPoints::add(verts[v].position, 0, 0.5, 0);
    }
#endif



}



bool OutputControllerEdgeFlipper::VertexRemovable(int vi) const {
    // a vertex is removable if it and all it's neighbors are NOT on any fronts

    const EdgeFlipperVertex &v = verts[vi];
    if (v.onfront)
	return false;

    for (std::set<int>::const_iterator iv=v.incidentVerticesLocal.begin();
	 iv != v.incidentVerticesLocal.end();
	 ++iv) {
	if (verts[*iv].onfront)
	    return false;
    }

    return true;
}

bool OutputControllerEdgeFlipper::TriangleFlippable(int ti) const {

    for (int i=0; i<3; i++) {
	if (verts[tris[ti].vi[i]].onfront)
	    return false;
    }

    return true;
}




void OutputControllerEdgeFlipper::FinishVertex(int v) {


    // free each triangle that is still here
    vector<int> ft;
    for (std::set<int>::iterator ti=verts[v].incidentTriangles.begin();
	 ti!=verts[v].incidentTriangles.end();
	 ++ti) {

	if (*ti >= 0)
	    ft.push_back(*ti);
    }

    // finish the 1-ring triangles
    for (unsigned ti=0; ti<ft.size(); ti++) {
	FinishTriangle(ft[ti]);
    }


    // remove v from it's neighbors as incident
    for (std::set<int>::iterator vi=verts[v].incidentVerticesLocal.begin();
	 vi!=verts[v].incidentVerticesLocal.end();
	 ++vi) {

	std::set<int>::iterator vi2 = verts[*vi].incidentVerticesLocal.find(v);
	if (vi2 != verts[*vi].incidentVerticesLocal.end())
	    verts[*vi].incidentVerticesLocal.erase(vi2);
    }



    for (std::set<int>::iterator ti=verts[v].incidentTriangles.begin();
	 ti!=verts[v].incidentTriangles.end();
	 ++ti) {

	if (*ti >= 0)
	    cerr<<"erasing triangle with incident verts??"<<endl;

    }

    // finish down stream
    child->FinalizeVertex(verts[v].index);


    global_to_local_vi.erase(global_to_local_vi.find(verts[v].index));
    FreeVert(v);

}



void OutputControllerEdgeFlipper::FinishTriangle(int t) {

    EdgeFlipperTriangle &tri = tris[t];

    // "remove" the triangle from its 3 vertices
    for (int i=0; i<3; i++) {
	EdgeFlipperVertex &v = verts[tri.vi[i]];
	std::set<int>::iterator ti = v.incidentTriangles.find(t);
	if (ti == v.incidentTriangles.end()) {
	    cerr<<"triangle not found in it's vertex!"<<endl;
	} else {
	    v.incidentTriangles.erase(ti);
	    v.incidentTriangles.insert(-t-1);
	}
    }


    // sometimes the we screw things up somehow
    if (verts[tri.vi[0]].index != verts[tri.vi[1]].index && 
	verts[tri.vi[1]].index != verts[tri.vi[2]].index &&
	verts[tri.vi[2]].index != verts[tri.vi[0]].index)

	child->AddTriangle(tri.index,
			   verts[tri.vi[0]].index,
			   verts[tri.vi[1]].index,
			   verts[tri.vi[2]].index);

    FreeTri(t);
}





int OutputControllerEdgeFlipper::AllocVert() {
    if (!available_verts.size()) {
	available_verts.push_back(verts.size());
	verts.resize(verts.size()+1);
    }

    int ret = available_verts.back();
    available_verts.resize(available_verts.size()-1);
    return ret;
}

int OutputControllerEdgeFlipper::AllocTri() {
    if (!available_tris.size()) {
	available_tris.push_back(tris.size());
	tris.resize(tris.size()+1);
    }

    int ret = available_tris.back();
    available_tris.resize(available_tris.size()-1);
    return ret;
}


void OutputControllerEdgeFlipper::FreeVert(int v) {
    available_verts.push_back(v);
}

void OutputControllerEdgeFlipper::FreeTri(int t) {
    available_tris.push_back(t);
}
