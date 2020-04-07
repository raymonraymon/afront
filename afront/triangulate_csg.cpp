
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
#include "triangulate_mesh.h"
#include "triangulate_mls.h"
#include "triangulate_csg.h"

#include "lls_wrapper.h"


MeshCSGGuidanceField::MeshCSGGuidanceField(int curv_sub, const TriangleMesh &mesh1, const TriangleMesh &mesh2, const vector<int> pointsides[2], vector< vector<Point3> > &curves, real_type rho, real_type min_step, real_type max_step, real_type reduction)
    : GuidanceField(rho, min_step, max_step, reduction), kdGetPoint() {


    // use the regular guidance field to find the point curvatures
    for (int i=0; i<2; i++) {
	const TriangleMesh &mesh = (i==0)?mesh1:mesh2;
	MeshGuidanceField gf(curv_sub, mesh, rho, min_step, max_step, reduction);

	for (unsigned j=0; j<mesh.verts.size(); j++) {
	    if (pointsides[i][j] > 0) {
		ideal_length.push_back(gf.IdealLength(j));
		kdGetPoint.allpoints.push_back(mesh.verts[j].point);
	    }
	}
    }


    // compute the curvature of the intersection curves
    for (unsigned c=0; c<curves.size(); c++) {

	for (unsigned cv=0; cv<curves[c].size(); cv++) {

	    vector<Point3> cpts;
	    vector<real_type> weights;

	    const int num=4;
	    for (int i=0; i<num*2+1; i++) {
		cpts.push_back(curves[c][(cv + curves[c].size() + i - num)%curves[c].size()]);
		weights.push_back(1);
	    }

	    ideal_length.push_back(MaxCurvatureToIdeal(edge_curvature(cpts, weights, cpts.size()/2)));
	    kdGetPoint.allpoints.push_back(curves[c][cv]);
	}

    }


    // setup the kdtree
    kdtree = new kdtree_type(10, Box3::make_union(mesh1.bounding_box(),mesh2.bounding_box()), kdGetPoint);
    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++)
        kdtree->Insert(i);
    kdtree->MakeTree();

}


MeshCSGGuidanceField::~MeshCSGGuidanceField() {
    delete kdtree;
}



void* MeshCSGGuidanceField::OrderedPointTraverseStart(const Point3 &p) {
	return new kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
}


bool MeshCSGGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = kdGetPoint(n);
    ideal = ideal_length[n];
    return true;
}

void MeshCSGGuidanceField::OrderedPointTraverseEnd(void *ctx) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    delete kdOrderedTraverse;
}



void MeshCSGGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad) {
    pts.resize(0);
    norms.resize(0);
    rad.resize(0);
}







Vector3 TriPointNormal(const Triangle3 &tri, const Point3 &p, const Vector3 &n0, const Vector3 &n1, const Vector3 &n2) {

    real_type b[3];
    tri.get_barycentric_coordinates(p, b[0], b[1], b[2]);

    Vector3 ret = b[0]*n0 + b[1]*n1 + b[2]*n2;
    ret.normalize();
    return ret;
}


void MeshCSGGuidanceField::GetIntersectionLoops(const TriangleMesh &m1, const TriangleMesh &m2, vector< vector<Point3> > &loops, vector< vector<Vector3> > &n1, vector< vector<Vector3> > &n2, vector<int> pointsides[2]) {

    vector<bool> vadded;
    gtb::tsurfel_set<real_type> startpoints;
    vector<Point3> endpoints;
    vector<Vector3> normals1, normals2;	// go with the startpoints, not the endpoints

    vector<int> kdfi;
    for (unsigned i=0; i<m1.faces.size(); i++) {
	kdfi.push_back(i);
    }

    TriangleMeshFaceTree kdtreebbox(m1);
    gtb::BoxKDTree<int, TriangleMeshFaceTree> kdtree(kdfi, kdtreebbox);

    vector< vector< std::pair<int,int> > > iverts[2];


    for (unsigned fi=0; fi<m2.faces.size(); fi++) {

	Box3 fbbox = Box3::bounding_box(m2.verts[m2.faces[fi].verts[0]].point,
					m2.verts[m2.faces[fi].verts[1]].point,
					m2.verts[m2.faces[fi].verts[2]].point);


	vector<int> possible;
	kdtree.GetIntersectedBoxes(kdtreebbox, fbbox, possible);
        
	for (int p=0; p<(int)possible.size(); p++) {


	    Triangle3 tri1(m1.verts[m1.faces[possible[p]].verts[0]].point, 
			   m1.verts[m1.faces[possible[p]].verts[1]].point, 
			   m1.verts[m1.faces[possible[p]].verts[2]].point);
	    Triangle3 tri2(m2.verts[m2.faces[fi].verts[0]].point, 
			   m2.verts[m2.faces[fi].verts[1]].point, 
			   m2.verts[m2.faces[fi].verts[2]].point);


	    Point3 int1, int2;
	    if (Triangle3::Intersection(tri1, tri2, int1, int2)) {
		startpoints.insert_vertex(int1);
		endpoints.push_back(int2);
		normals1.push_back(TriPointNormal(tri1, int1, m1.verts[m1.faces[possible[p]].verts[0]].normal,
						  m1.verts[m1.faces[possible[p]].verts[1]].normal,
						  m1.verts[m1.faces[possible[p]].verts[2]].normal));
		normals2.push_back(TriPointNormal(tri1, int1, m2.verts[m2.faces[fi].verts[0]].normal,
						  m2.verts[m2.faces[fi].verts[1]].normal,
						  m2.verts[m2.faces[fi].verts[2]].normal));
		vadded.push_back(false);

		vector< std::pair<int,int> > tiverts[2];

		for (int i=0; i<3; i++) {

		    real_type t;
		    if (tri1.intersect_segment(tri2[i], tri2[(i+1)%3], t)) {
			if (tri1.normal().dot(tri2[i]) > tri1.normal().dot(tri1[0])) {
			    tiverts[1].push_back(std::pair<int,int>(m2.faces[fi].verts[i], 1));
			    tiverts[1].push_back(std::pair<int,int>(m2.faces[fi].verts[(i+1)%3], -1));
			} else {
			    tiverts[1].push_back(std::pair<int,int>(m2.faces[fi].verts[i], -1));
			    tiverts[1].push_back(std::pair<int,int>(m2.faces[fi].verts[(i+1)%3], 1));
			}
		    }

		    if (tri2.intersect_segment(tri1[i], tri1[(i+1)%3], t)) {
			if (tri2.normal().dot(tri1[i]) > tri2.normal().dot(tri2[0])) {
			    tiverts[0].push_back(std::pair<int,int>(m1.faces[possible[p]].verts[i], 1));
			    tiverts[0].push_back(std::pair<int,int>(m1.faces[possible[p]].verts[(i+1)%3], -1));
			} else {
			    tiverts[0].push_back(std::pair<int,int>(m1.faces[possible[p]].verts[i], -1));
			    tiverts[0].push_back(std::pair<int,int>(m1.faces[possible[p]].verts[(i+1)%3], 1));
			}
		    }

		}

		iverts[0].push_back(tiverts[0]);
		iverts[1].push_back(tiverts[1]);

	    }
	}

    }


    gtb::ss_kdtree<surfel_set> kd(startpoints);
    bool useall=false;
    while (1) {

	// find a vert that hasn't been used
	int starti=0;
	for (starti=0; starti<(int)vadded.size(); starti++) {
	    if (!vadded[starti]) break;
	}

	// all verts were added
	if (starti>=(int)vadded.size())	break;

	vector<int> iloop;
	iloop.push_back(starti);

	vadded[starti] = true;
	int last = starti;

	while (true) {

            int next = kd.tree->FindMin(endpoints[last]);
	    if (next == starti) break;

	    iloop.push_back(next);

	    vadded[next] = true;
	    last = next;

	}


	bool use = false;

	// use interface to decide which intersection curves to use
	if (!useall) {
	    dbgClear();
	    for (unsigned i=0; i<iloop.size(); i++) {
		DbgPoints::add(startpoints.vertex(iloop[i]), 1, 0, 0);
	    }

	    while (1) {
		cerr<<"use this loop [y/n/a]?"<<endl;
		redrawAndWait();

		//				int lastkey = GetUILastKey();
		int lastkey = 'a';
		if (lastkey=='n' || lastkey=='N') {
		    use = false; break;
		}

		if (lastkey=='y' || lastkey=='Y') {
		    use = true; break;
		}

		if (lastkey=='a' || lastkey=='A') {
		    useall = true; break;
		}

	    }

	}

	if (useall || use) {
	    vector<Point3> ploop;
	    vector<Vector3> n1loop;
	    vector<Vector3> n2loop;
	    for (unsigned i=0; i<iloop.size(); i++) {
		ploop.push_back(startpoints.vertex(iloop[i]));
		n1loop.push_back(normals1[iloop[i]]);
		n2loop.push_back(normals2[iloop[i]]);

		for (unsigned j=0; j<iverts[0][iloop[i]].size(); j++) {
		    pointsides[0][iverts[0][iloop[i]][j].first] = iverts[0][iloop[i]][j].second;
		}
		for (unsigned j=0; j<iverts[1][iloop[i]].size(); j++) {
		    pointsides[1][iverts[1][iloop[i]][j].first] = iverts[1][iloop[i]][j].second;
		}
	    }
	    loops.push_back(ploop);
	    n1.push_back(n1loop);
	    n2.push_back(n2loop);
	}
    }
}



void MeshCSGGuidanceField::FloodPointSides(const TriangleMesh &m, vector<int> &sides) {

    vector<int> toprocess;
    for (unsigned i=0; i<m.verts.size(); i++) {
	if (sides[i] != 0)
	    toprocess.push_back(i);
    }

    while (toprocess.size()) {

	int from = toprocess.back();
	toprocess.pop_back();

	for (TriangleMesh::VertexVertexIteratorI vi(m, from); !vi.done(); ++vi) {

	    int nbr = *vi;
	    if (sides[nbr] == 0) {
		sides[nbr] = sides[from];
		toprocess.push_back(nbr);
	    }
	}
    }
}


void MeshCSGGuidanceField::TrimPointSides(const TriangleMesh &mesh, vector<int> &sides) {

    real_type max_dist = mesh.bounding_box().diagonal_length() * 0.06;

    gtb::fast_pq< std::pair<real_type, int> > pq;

    vector<real_type> dists(sides.size());
    for (unsigned i=0; i<dists.size(); i++) {
	dists[i] = -1;

	for (TriangleMesh::VertexVertexIteratorI vi(mesh, i); !vi.done(); ++vi) {

	    if (sides[*vi] != sides[i]) {
		pq.push(std::pair<real_type,int>(0,i));
	    }
	}
    }


    while (!pq.empty()) {

	if (-pq.top().first>max_dist) break;

	int v = pq.top().second;

	if (dists[v] >= 0) {
	    pq.pop();
	    continue;	// already set
	}

	if (-pq.top().first > 3*MaxStepLength(mesh.verts[v].point) / (1-reduction)) {
	    dists[v] = 0;
	    pq.pop();
	    continue;
	}

	dists[v] = -pq.top().first;
	sides[v] *= 2;

	pq.pop();


	for (TriangleMesh::VertexVertexIteratorI vi(mesh, v); !vi.done(); ++vi) {

	    if (dists[*vi] < 0) {
		pq.push(std::pair<real_type,int>(-(dists[v] + Point3::distance(mesh.verts[v].point, mesh.verts[*vi].point)),*vi));
	    }
	}
    }

    FixPointSides(mesh, sides);
}




void MeshCSGGuidanceField::FixPointSides(const TriangleMesh &mesh, vector<int> &sides) {


    std::set<int> tofix;

    for (unsigned i=0; i<mesh.verts.size(); i++) {

	if (sides[i] != 1) continue;

	for (TriangleMesh::VertexVertexIteratorI vi(mesh, i); !vi.done(); ++vi) {

	    if (sides[*vi] != 1) {
		tofix.insert(i);
		break;
	    }
	}
    }


    while (!tofix.empty()) {

	int v = *(tofix.begin());
	tofix.erase(tofix.begin());

	if (sides[v] != 1) continue;	// already moved

	// see if it really needs fixing
	static_vector(int, ring);
	bool edge = mesh.Vert1Ring(v, ring);

	int stepouts=0;
	int numout=0;
	for (unsigned i=0; i<(edge?ring.size()-1:ring.size()); i++) {
	    if (sides[ring[i]]==1 && sides[ring[(i+1)%ring.size()]]!=1)
		stepouts++;
	    if (sides[ring[i]]==1)
		numout++;
	}

	// move it into the 2 area
	if (stepouts!=1 || numout==1) {
	    cerr<<"fixing vert side"<<endl;


	    if (0) {
		dbgClear();
		for (unsigned i=0; i<mesh.verts.size(); i++) {
		    if (sides[i] == -1) {
			DbgPoints::add(mesh.verts[i].point, 0.5, 0, 0);
		    } else if (sides[i] == -2) {
			DbgPoints::add(mesh.verts[i].point, 1, 0, 0);
		    } else if (sides[i] == 1) {
			DbgPoints::add(mesh.verts[i].point, 0, 0.5, 0);
		    } else if (sides[i] == 2) {
			DbgPoints::add(mesh.verts[i].point, 0, 1, 0);
		    }
		}
		DbgPoints::add(mesh.verts[v].point, 0, 0, 1);
		redrawAndWait(' ');
	    }



	    sides[v] = 2;

	    for (unsigned i=0; i<ring.size(); i++) {
		if (sides[ring[i]] == 1)
		    tofix.insert(ring[i]);
	    }
	}
    }
}


void MeshCSGGuidanceField::FrontsFromSides(const TriangleMesh &mesh, vector<int> &sides, vector< vector<Point3> > &fpoints, vector< vector<Vector3> > &fnormals) {

    std::set<int> onfront;
    for (unsigned i=0; i<mesh.verts.size(); i++) {
	if (sides[i] != 1) continue;
	for (TriangleMesh::VertexVertexIteratorI vi(mesh, i); !vi.done(); ++vi) {
	    if (sides[*vi] != 1) {
		onfront.insert(i);
		break;
	    }
	}
    }

    while (!onfront.empty()) {

	vector<int> front;

	front.push_back(*onfront.begin());
	onfront.erase(onfront.begin());

	while (1) {

	    static_vector(int, ring);
	    bool edge = mesh.Vert1Ring(front.back(), ring);

	    int toadd=-1;
	    for (unsigned i=0; i<(edge?ring.size()-1:ring.size()); i++) {
		if (sides[ring[i]]==1 && sides[ring[(i+1)%ring.size()]]!=1) {
		    if (toadd!=-1) {
			cerr<<"more than one stepout!"<<endl;
			BREAK;
		    }
		    toadd = ring[i];
		}
	    }

	    if (toadd == -1) {
		cerr<<"didn't know where to go!"<<endl;
		BREAK;
	    }

	    if (toadd == front.front()) break;

	    onfront.erase(onfront.find(toadd));
	    front.push_back(toadd);
	}

	fpoints.push_back(vector<Point3>());
	fnormals.push_back(vector<Vector3>());

	for (unsigned i=0; i<front.size(); i++) {
	    fpoints.back().push_back(mesh.verts[front[i]].point);
	    fnormals.back().push_back(mesh.verts[front[i]].normal);
	}
    }
}




// Floater mean-value weights
real_type floaterMV(const vector<Vector3> &vnei, vector<real_type> &w) {
    unsigned n = vnei.size();
    // flatten neighborhood using Z^alpha map; keep results in polar coordinates (r,phi)
    vector<real_type> r,phi,ta_2; r.resize(n); phi.resize(n); ta_2.resize(n);
    real_type sum_phi=0;
    for (unsigned i=0; i<n; i++){
        r[i] = vnei[i].length();
        if(i) sum_phi += acossafe( vnei[i].dot(vnei[i-1]) / (r[i] * r[i-1]));
        phi[i] = sum_phi;
    }
    sum_phi += acossafe( vnei[0].dot(vnei[n-1]) / (r[0] * r[n-1]));
    sum_phi = 2 * (real_type)M_PI / sum_phi;
    for (unsigned i=0; i<n; i++) { phi[i] *= sum_phi; } // normalize angles in order to bring sum to 2PI
    for (unsigned i=0; i<n; i++) {
	real_type alpha=phi[(i+1)%n]-phi[i];
	ta_2[i] = tan(alpha*0.5);
    }  // find angles between the spokes
    //the original Z^alpha map also scaled r's using an exponent dependent on sum_phi
    w.resize(n); real_type sumw=0;
    for (unsigned i=0; i<n; i++) { w[i] = (ta_2[(i+n-1)%n] + ta_2[i]) / r[i]; sumw+=w[i]; }
    //caller function still needs to do normalization ...
    return sumw;
}


void MeshCSGGuidanceField::BlendGuidance(bool blend, const TriangleMesh* meshes[2], vector<int> psides[2]) {

    if (!blend) {
	return;
    }

    int pindex=0;

    for (int m=0; m<2; m++) {

	const TriangleMesh &mesh = *meshes[m];
	vector<int> sides = psides[m];

	// count the number of verts we will be solving for
	// side=1: what we want
	// side=2: the existing points
	// side<0: the deleted points

	vector<int> mesh_to_lls(mesh.verts.size()), lls_to_mesh;

	vector<real_type> priorities(mesh.verts.size());
	for (unsigned i=0; i<sides.size(); i++) {
	    //			cerr<<sides[i]<<endl;
	    if (sides[i] == 2) {
		mesh_to_lls[i] = lls_to_mesh.size();
		lls_to_mesh.push_back(i);
		priorities[i] = -1;
	    } else if (sides[i]==1) {
		mesh_to_lls[i]=-1;
		priorities[i] = 1;
	    } else {
		mesh_to_lls[i]=-1;
		priorities[i] = 0;
	    }
	}


	LLSWrapper<real_type> lls(lls_to_mesh.size(), lls_to_mesh.size());

	vector<int> nbrs;
	vector<Vector3> vnei;
	vector<real_type> weights;
	for (unsigned r=0; r<lls_to_mesh.size(); r++) {

	    nbrs.clear();
	    for (TriangleMesh::VertexVertexIteratorI vi(mesh, lls_to_mesh[r]); !vi.done(); ++vi) {
		nbrs.push_back(*vi);
	    }

	    // check if we're on a boundary
	    bool onbound=false;
	    for (unsigned i=0; i<nbrs.size(); i++) {
		if (sides[nbrs[i]] < 1) {
		    lls.InsertA(r,r, 1);
		    lls.InsertB(r,0);
		    lls.InsertX(r,0);
		    onbound=true; break;
		}
		if (sides[nbrs[i]] == 1) {
		    lls.InsertA(r,r, 1);
		    lls.InsertB(r,1);
		    lls.InsertX(r,1);
		    onbound=true; break;
		}
	    }

	    if (onbound) continue;


	    // use mean value weights
	    vnei.clear();
	    for (unsigned i=0; i<nbrs.size(); i++) {
		vnei.push_back(mesh.verts[nbrs[i]].point - mesh.verts[lls_to_mesh[r]].point);
	    }

	    real_type wsumi = 1 / floaterMV(vnei, weights);

	    lls.InsertA(r,r,-1);
	    for (unsigned i=0; i<nbrs.size(); i++) {
		lls.InsertA(r,mesh_to_lls[nbrs[i]], wsumi*weights[i]);
	    }

	}


	lls.Solve();

	cerr<<"lls size: "<<lls_to_mesh.size()<<endl;
	cerr<<"meshsize: "<<mesh.verts.size()<<endl;



	for (unsigned r=0; r<lls_to_mesh.size(); r++) {
	    priorities[lls_to_mesh[r]] = lls.GetX(r);
	}

#if 1
	{
	    dbgClear();
	    for (unsigned f=0; f<mesh.faces.size(); f++) {
		DbgPoly::add(mesh.verts[mesh.faces[f].verts[0]].point, Point3(priorities[mesh.faces[f].verts[0]], 0, 1-priorities[mesh.faces[f].verts[0]]),
			     mesh.verts[mesh.faces[f].verts[1]].point, Point3(priorities[mesh.faces[f].verts[1]], 0, 1-priorities[mesh.faces[f].verts[1]]),
			     mesh.verts[mesh.faces[f].verts[2]].point, Point3(priorities[mesh.faces[f].verts[2]], 0, 1-priorities[mesh.faces[f].verts[2]]));

	    }
	    redrawAndWait(' ');
	}
#endif


	for (unsigned i=0; i<sides.size(); i++) {

	    if (sides[i] > 0) {

		if (mesh_to_lls[i] >= 0) {

		    // blend the curvature radius with the existing edge lengths

		    real_type ave_elen=0;
		    int nnbrs=0;

		    for (TriangleMesh::VertexVertexIteratorI vi(mesh, i); !vi.done(); ++vi) {
			ave_elen += Point3::distance(mesh.verts[*vi].point, mesh.verts[i].point);
			nnbrs++;
		    }
		    ave_elen /= nnbrs;

		    ideal_length[pindex] = priorities[i]*ave_elen + (1-priorities[i])*ideal_length[pindex];
		} else {
		    ideal_length[pindex] = 100000;
		}

		pindex++;

	    }
	}
    }

    cerr<<"after blend: "<<pindex<<endl;

}





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////



MeshPSCSGGuidanceField::MeshPSCSGGuidanceField(int curv_sub, const TriangleMesh &mesh, const surfel_set &points, CProjection &psprojector, vector< vector<Point3> > &curves, real_type rho, real_type min_step, real_type max_step, real_type reduction, int adamson)
    : GuidanceField(rho, min_step, max_step, reduction), kdGetPoint() {



    // use the regular guidance field to find the point curvatures
    {
	MeshGuidanceField gf(curv_sub, mesh, rho, min_step, max_step, reduction);

	for (unsigned j=0; j<mesh.verts.size(); j++) {
	    if (true) { //pointsides[0][j] > 0) {
		ideal_length.push_back(gf.IdealLength(j));
		kdGetPoint.allpoints.push_back(mesh.verts[j].point);
	    }
	}
    }

    {
	SmoothMLSGuidanceField gf(psprojector, rho, min_step, max_step, reduction, adamson);

	for (unsigned j=0; j<points.size(); j++) {
	    if (true) { //pointsides[1][j] > 0) {
		ideal_length.push_back(gf.IdealLength(j));
		kdGetPoint.allpoints.push_back(mesh.verts[j].point);
	    }
	}
    }



    // compute the curvature of the intersection curves
    for (unsigned c=0; c<curves.size(); c++) {
	for (unsigned cv=0; cv<curves[c].size(); cv++) {

	    vector<Point3> cpts;
	    vector<real_type> weights;

	    const int num=4;
	    for (int i=0; i<num*2+1; i++) {
		cpts.push_back(curves[c][(cv + curves[c].size() + i - num)%curves[c].size()]);
		weights.push_back(1);
	    }

	    ideal_length.push_back(MaxCurvatureToIdeal(edge_curvature(cpts, weights, cpts.size()/2)));
	    kdGetPoint.allpoints.push_back(curves[c][cv]);
	}
    }


    // setup the kdtree
    kdtree = new kdtree_type(10, Box3::make_union(mesh.bounding_box(), points.bounding_box()), kdGetPoint);
    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++)
        kdtree->Insert(i);
    kdtree->MakeTree();

}


MeshPSCSGGuidanceField::~MeshPSCSGGuidanceField() {
    delete kdtree;
}


void* MeshPSCSGGuidanceField::OrderedPointTraverseStart(const Point3 &p) {
    return new kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
}


int MeshPSCSGGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    if (kdOrderedTraverse->empty()) return false;
    return kdOrderedTraverse->GetNext(squared_dist);
}

bool MeshPSCSGGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = kdGetPoint(n);
    ideal = ideal_length[n];
    return true;
}


void MeshPSCSGGuidanceField::OrderedPointTraverseEnd(void *ctx) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    delete kdOrderedTraverse;
}


void MeshPSCSGGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad) {
    pts.resize(0);
    norms.resize(0);
    rad.resize(0);
}


void MeshPSCSGGuidanceField::ProjectPointOntoIntersection(MeshProjector &mp, SmoothMLSProjector &psp,
							  const Point3 &from, const Point3 *last, real_type dist, Point3 &to, Vector3 &mnorm, Vector3 &pnorm) {

    to = from;

    int iterations = 0;

    while (1) {
	iterations++;

	// project onto the mesh
	Point3 mnext;
	mp.ProjectPoint(to, mnext, mnorm);
	if (last) {
	    Vector3 dir = mnext - *last;
	    dir.normalize();
	    mnext = *last + dist*dir;
	}


#if 0
	{
	    dbgClear();
	    if (last) DbgPoints::add(*last, 0, 1, 0);
	    DbgPoints::add(mnext, 1, 0, 0);
	    redrawAndWait(' ');
	}
#endif


	// project onto the pointset
	Point3 pnext;
	psp.ProjectPoint(mnext, pnext, pnorm);
	if (last) {
	    Vector3 dir = pnext - *last;
	    dir.normalize();
	    mnext = *last + dist*dir;
	}


#if 0
	{
	    dbgClear();
	    if (last) DbgPoints::add(*last, 0, 1, 0);
	    DbgPoints::add(pnext, 1, 0, 0);
	    redrawAndWait(' ');
	}
#endif


	if (iterations > 10)
	    break;
	if (Point3::distance(to, pnext) < 0.000001)
	    break;

	to = pnext;
    }
}


void MeshPSCSGGuidanceField::GetIntersectionLoops(MeshProjector &mp, SmoothMLSProjector &psp,
						  vector< vector<Point3> > &loops, vector< vector<Vector3> > &mnorms, vector< vector<Vector3> > &pnorms) {

    Point3 p;

    // big ugly hack so I don't have to keep selecting the points over and over
    Point3 selpts[] = {
	Point3(-0.21770936250686646,0.0017753387801349163,-0.025098925456404686),
	Point3(-0.19265307486057281,0.10258496552705765,-0.025064017623662949),
	Point3(-0.17150920629501343,0.13518926501274109,-0.0250661950558424),
	Point3(-0.13799993693828583,0.1679045557975769,-0.024924393743276596),
	Point3(-0.10503314435482025,0.19052994251251221,-0.024971576407551765),
	Point3(0.00086189928697422147,0.21736522018909454,-0.024933561682701111),
	Point3(0.10065957903862,0.19173747301101685,-0.024877807125449181),
	Point3(0.1387404203414917,0.1683715432882309,-0.024997999891638756),
	Point3(0.166483074426651,0.13879218697547913,-0.024914190173149109),
	Point3(0.19156813621520996,0.10502850264310837,-0.025036506354808807),
	Point3(0.21675296127796173,-0.0003444627218414098,-0.025008678436279297),
	Point3(0.19217985868453979,-0.099696114659309387,-0.025007132440805435),
	Point3(0.17041212320327759,-0.13730905950069427,-0.02548324316740036),
	Point3(0.14048014581203461,-0.16731677949428558,-0.025159848853945732),
	Point3(0.10266976058483124,-0.19194731116294861,-0.025125492364168167),
	Point3(0.00062523921951651573,-0.21382820606231689,-0.024954486638307571),
	Point3(-0.10081255435943604,-0.18937049806118011,-0.024979390203952789),
	Point3(-0.13664749264717102,-0.16828446090221405,-0.025129873305559158),
	Point3(-0.16698750853538513,-0.13863109052181244,-0.02513224259018898),
	Point3(-0.19095830619335175,-0.10008219629526138,-0.024988966062664986),
	Point3(-0.27817359566688538,-0.0016937416512519121,-0.03019353374838829),
	Point3(-0.24604599177837372,0.12733061611652374,-0.029946429654955864),
	Point3(-0.21648073196411133,0.17285509407520294,-0.02994309738278389),
	Point3(-0.17454756796360016,0.21677061915397644,-0.030023155733942986),
	Point3(-0.12787738442420959,0.24664844572544098,-0.029966993257403374),
	Point3(-0.0031186200212687254,0.27720457315444946,-0.029851516708731651),
	Point3(0.12436946481466293,0.24750843644142151,-0.0298641137778759),
	Point3(0.17439842224121094,0.21549306809902191,-0.029891205951571465),
	Point3(0.2163376659154892,0.17340311408042908,-0.029925864189863205),
	Point3(0.2471354752779007,0.12469444423913956,-0.02992076613008976),
	Point3(0.27787220478057861,0.0071234251372516155,-0.030091237276792526),
	Point3(0.24729768931865692,-0.12573428452014923,-0.030009320005774498),
	Point3(0.21625505387783051,-0.1730244904756546,-0.030011618509888649),
	Point3(0.17580631375312805,-0.21388167142868042,-0.030024290084838867),
	Point3(0.12717312574386597,-0.24709805846214294,-0.030163409188389778),
	Point3(0.0073157469742000103,-0.27722537517547607,-0.030166942626237869),
	Point3(-0.12562718987464905,-0.24776168167591095,-0.030323376879096031),
	Point3(-0.17502468824386597,-0.21622093021869659,-0.030331624671816826),
	Point3(-0.2188437283039093,-0.17430791258811951,-0.030237002298235893),
	Point3(-0.2475532591342926,-0.12683656811714172,-0.030277537181973457),
    };


    real_type step = 0.001;

    vector<Point3> loop;
    vector<Vector3> tmnorms, tpnorms;

    int sp=0;

    while (1) {

#if 0 	// either take points from the gui, or from the array above
	cerr<<"select point to start projection from"<<endl;

	bool done=false;
	while (1) {
	    redrawAndWait(' ');
	    int key = GetUILastKey();
	    if (key == ' ') {
		if (getUISelectPoint(p, win)) {
		    break;
		}
	    } else if (key == 'b' || key == 'B') {
		done=true; break;
	    }
	}
	if (done) break;

	cerr<<p<<endl;
	continue;
#else
	if (sp >= sizeof(selpts)/sizeof(selpts[0])) break;
	p = selpts[sp]; sp++;
#endif



#if 0
	{
	    dbgClear();
	    DbgPoints::add(p, 1, 0, 0);
	    redrawAndWait(' ');
	}
#endif

	loop.clear();
	tmnorms.clear();
	tpnorms.clear();

	
	Point3 to;
	Vector3 mnorm, pnorm;


	ProjectPointOntoIntersection(mp, psp, p, NULL, 0, to, mnorm, pnorm);
	loop.push_back(to);
	tmnorms.push_back(mnorm);
	tpnorms.push_back(pnorm);


	while (1) {

	    Vector3 dir = tmnorms.back().cross(tpnorms.back());
	    dir.normalize();
	    Point3 from = loop.back() + step * dir;

	    {
		dbgClear();
				
		{
		    vector<Point3> line(2);
		    line[0] = loop.back();
		    line[1] = line[0] + (real_type)0.02 * tmnorms.back();
		    DbgPLines::add(line, 2);
		    line[1] = line[0] + (real_type)0.02 * tpnorms.back();
		    DbgPLines::add(line, 1);
		}

		for (unsigned i=0; i<loop.size(); i++) {
		    DbgPoints::add(loop[i], 1, 0, 0);
		}
		DbgPoints::add(from, 0, 1, 0);
		redrawAndWait();
	    }


	    ProjectPointOntoIntersection(mp, psp, from, &loop.back(), step, to, mnorm, pnorm);
	    loop.push_back(to);
	    tmnorms.push_back(mnorm);
	    tpnorms.push_back(pnorm);

	    {
		dbgClear();
		for (unsigned i=0; i<loop.size(); i++) {
		    DbgPoints::add(loop[i], 1, 0, 0);
		}
		redrawAndWait();
	    }


	    if (loop.size()>5 && Point3::distance(loop.front(), loop.back())<2*step)
		break;
	}

	loops.push_back(loop);
	mnorms.push_back(tmnorms);
	pnorms.push_back(tpnorms);
    }
}




void MeshPSCSGGuidanceField::InitPointSides(const TriangleMesh &m, const surfel_set &points, vector<int> &sides, const vector< vector<Point3> > &loops) {


    sides.resize(m.verts.size());
    for (int i=0; i<(int)m.verts.size(); i++) {
        sides[i] = 0;
    }


    // get all the mesh points to check
    for (int l=0; l<(int)loops.size(); l++) {
	for (int ll=0; ll<(int)loops[l].size(); ll++) {
	    void *ctx = OrderedPointTraverseStart(loops[l][ll]);
	    while(1) {
		real_type sdist;
		int next = OrderedPointTraverseNext(ctx, sdist);
		if (next < 0) break;
		if (sqrt(sdist) > 0.003) break;
		if (next < (int)m.verts.size()) {
		    sides[next] = 10;
		}
	    }
	    OrderedPointTraverseEnd(ctx);
	}
    }
    cerr<<"got mesh points to check"<<endl;

    dbgClear();
    for (int v=0; v<(int)m.verts.size(); v++) {
	if (sides[v]<0) DbgPoints::add(m.verts[v].point, 0,0,1);
	if (sides[v]>0) DbgPoints::add(m.verts[v].point, 1,0,0);
    }
    redrawAndWait(' ', true);
    dbgClear();


    // check them against the pointset points
    for (unsigned v=0; v<sides.size(); v++) {

	if (sides[v] != 10) continue;
	sides[v] = 0;

	void *ctx = OrderedPointTraverseStart(m.verts[v].point);

	while(1) {
	    real_type sdist;
	    int next = OrderedPointTraverseNext(ctx, sdist);
	    if (next < 0) break;
	    //			if (sqrt(sdist) > 0.01) break;
	    if (next >= (int)m.verts.size() && next <(int)m.verts.size()+(int)points.size()) {

		if ((m.verts[v].point - points.vertex(next-m.verts.size())).dot(points.normal(next-m.verts.size())) > 0)
		    sides[v] -= 1;
		else
		    sides[v] += 1;

		break;
	    }
	}
	OrderedPointTraverseEnd(ctx);

	if (sides[v] < 0) sides[v]=-2;
	else if (sides[v] > 0) sides[v]=2;
	else cerr<<"tie in point side voting!"<<endl;


	if (sides[v]==-2) DbgPoints::add(m.verts[v].point, 0,0,1);
	if (sides[v]== 2) DbgPoints::add(m.verts[v].point, 1,0,0);
    }

    cerr<<"done checking "<<endl;


    dbgClear();
    for (int v=0; v<(int)m.verts.size(); v++) {
	if (sides[v]==-2) DbgPoints::add(m.verts[v].point, 0,0,1);
	if (sides[v]== 2) DbgPoints::add(m.verts[v].point, 1,0,0);
    }
    redrawAndWait(' ', true);
}



// figure out which side of the torus the mesh vertices are on
void MeshPSCSGGuidanceField::HackPointSides(const TriangleMesh &m, vector<int> &sides) {
    cerr<<"doing HUUUUGGGGGEEE hack!!"<<endl;

    sides.resize(m.verts.size());
    for (int i=0; i<(int)m.verts.size(); i++) {
        sides[i] = 0;
    }


    for (int v=0; v<(int)m.verts.size(); v++) {
	
	Vector3 pv = m.verts[v].point - Point3(0,0,0);
	Vector3 zv=pv; zv[2]=0;
	zv.normalize();

	zv *= 0.25;	// torus radius

	real_type dist = (pv-zv).length();

	if (dist > 0.04) {
	    sides[v] = dist > 0.06 ? 1 : 2;
	} else {
	    sides[v] = dist < 0.02 ?  -1 :  -2;
	}
    }


    dbgClear();
    for (int v=0; v<(int)m.verts.size(); v++) {
	if (sides[v]==-2) DbgPoints::add(m.verts[v].point, 0,0,1);
	if (sides[v]== 2) DbgPoints::add(m.verts[v].point, 1,0,0);
	if (sides[v]==-1) DbgPoints::add(m.verts[v].point, 0,0,0.5);
	if (sides[v]== 1) DbgPoints::add(m.verts[v].point, 0.5,0,0);
    }
    redrawAndWait(' ', true);
}

