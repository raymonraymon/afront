
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
#include "triangulate_tet.h"

#include <lib/mlslib/NR/nr.h>	// for zbrent?!?!?


// define this to use the intersection of the ring and the level set as the projection
#define TET_CIRCULAR_PROJECTION
const real_type theta_step = M_PI_2/10;


extern int curvature_sub;
extern int eval_sub;


const real_type weight_scale = 0.2;
const int max_query_knn = 2500;
const int desired_knn = 250;
//const int desired_knn = 500;
const real_type cache_overestimate = 1.05;

const int nielson_gradient_smoothpasses = 3;

typedef long double projection_type;
typedef long double eval_type;
typedef long double normal_type;
typedef long double curvature_type;
typedef double curvature_comp_type;  // can't do long double since gtb::mat3 isn't smart enough


const bool lls_normal_equations = true;
//const gtb::poly_lls_solver lls_solver = gtb::POLY_LLS_SVD;
//const gtb::poly_lls_solver lls_solver = gtb::POLY_LLS_CHOLESKY;
const gtb::poly_lls_solver lls_solver = gtb::POLY_LLS_HOUSEHOLDER;

const real_type zbrent_tol = (real_type)1e-3;
const real_type step_dist_tol = (real_type)1e-5;
const real_type iso_dist_tol = (real_type)1e-5;


extern bool allow_outside;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// define this stuff so we only have to promote the type to Gradient<> when needed

template <typename T, typename P>
    gtb::tPoint3<T> gp(const gtb::tPoint3<P> &x) { return gtb::tPoint3<T>((T)x[0],(T)x[1],(T)x[2]); }

template <typename T, typename P>
    gtb::tVector3<T> gv(const gtb::tVector3<P> &x) { return gtb::tVector3<T>((T)x[0],(T)x[1],(T)x[2]); }


template <typename T, typename P>
    gtb::tPoint2<T> gp(const gtb::tPoint2<P> &x) {
    gtb::tPoint2<T> ret;
    ret[0] = x[0];
    ret[1] = x[1];
    return ret;
}

template <typename T, typename P>
    gtb::tVector2<T> gv(const gtb::tVector2<P> &x) {
    gtb::tVector2<T> ret;
    ret[0] = x[0];
    ret[1] = x[1];
    return ret;
}



template <class T>
    T cross(const gtb::tVector2<T> &v1, const gtb::tVector2<T> &v2)
{
    return (v1[0] * v2[1]) - (v1[1] * v2[0]);
}



template<typename T>
    gtb::tVector3<T> gsub1(const gtb::tPoint3<T> &lhs, const Point3 &rhs) {
    return gtb::tVector3<T>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2]);
}
template<typename T>
    gtb::tVector3<T> gsub2(const Point3 &lhs, const gtb::tPoint3<T> &rhs) {
    return gtb::tVector3<T>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2]);
}
template<typename T>
    gtb::tVector2<T> gsub1(const gtb::tPoint2<T> &lhs, const Point2 &rhs) {
    return gtb::tVector2<T>(lhs[0]-rhs[0], lhs[1]-rhs[1]);
}
template<typename T>
    gtb::tVector2<T> gsub2(const Point2 &lhs, const gtb::tPoint2<T> &rhs) {
    return gtb::tVector2<T>(lhs[0]-rhs[0], lhs[1]-rhs[1]);
}

template<typename T>
    gtb::tPoint2<T> gadd(const Point2 &lhs, const gtb::tVector2<T> &rhs) {
    return gtb::tPoint2<T>(lhs[0]+rhs[0], lhs[1]+rhs[1]);
}

template<typename T>
    void gadd_scaled(gtb::tPoint3<T> &lhs, const Point3 &rhs, const T &scale) {
    lhs[0] += rhs[0]*scale;	lhs[1] += rhs[1]*scale;	lhs[2] += rhs[2]*scale;
}

template<typename T>
    T gdot1(const gtb::tVector3<T> &lhs, const Vector3 &rhs) {
    return (lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2]);
}
template<typename T>
    T gdot2(const Vector3 &lhs, const gtb::tVector3<T> &rhs) {
    return (lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2]);
}
template<typename T>
    T gdot1(const gtb::tVector2<T> &lhs, const Vector2 &rhs) {
    return (lhs[0]*rhs[0] + lhs[1]*rhs[1]);
}
template<typename T>
    T gdot2(const Vector2 &lhs, const gtb::tVector2<T> &rhs) {
    return (lhs[0]*rhs[0] + lhs[1]*rhs[1]);
}


template<typename T>
    gtb::tVector3<T> gscale(const Vector3 &x, const T &scale) {
    return gtb::tVector3<T>(x[0]*scale, x[1]*scale, x[2]*scale);
}
template<typename T>
    gtb::tVector3<T> gscale(const T &scale, const Vector3 &x) {
    return gtb::tVector3<T>(x[0]*scale, x[1]*scale, x[2]*scale);
}
template<typename T>
    gtb::tVector2<T> gscale(const Vector2 &x, const T &scale) {
    return gtb::tVector2<T>(x[0]*scale, x[1]*scale);
}
template<typename T>
    gtb::tVector2<T> gscale(const T &scale, const Vector2 &x) {
    return gtb::tVector2<T>(x[0]*scale, x[1]*scale);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////




bool TetMesh::ReadOFF(const char *filename) {

    Clear();

    real_type min_value = 1e34;
    real_type max_value =-1e34;


    FILE *f = fopen(filename, "r");
    if (!f) return false;

    int nverts, ntets;
    if (fscanf(f, "%d %d\n", &nverts, &ntets) != 2) { fclose(f); Clear(); return false; }

    verts.resize(nverts);
    tets.resize(ntets);

    for (int i=0; i<nverts; i++) {
	float v[4] = { 0, 0, 0, 0};
	int r = fscanf(f, "%f %f %f %f\n", &v[0],&v[1],&v[2],&v[3]);
	if (r !=4 && r!=3) { fclose(f); Clear(); return false; }
	verts[i] = Vert(Point3(v[0],v[1],v[2]), v[3]);
	min_value = std::min(min_value, (real_type)v[3]);
	max_value = std::max(max_value, (real_type)v[3]);
    }
    

    for (int i=0; i<ntets; i++) {
	if (fscanf(f, "%d %d %d %d\n",
		   &tets[i].vi[0], &tets[i].vi[1], &tets[i].vi[2], &tets[i].vi[3]) != 4) { fclose(f); Clear(); return false; }

	// we want a consistent ordering of the tets
	real_type area = (verts[tets[i].vi[2]].point-verts[tets[i].vi[0]].point).cross(
										       verts[tets[i].vi[1]].point-verts[tets[i].vi[0]].point).dot(
																		  verts[tets[i].vi[3]].point-verts[tets[i].vi[0]].point);
	if (area < 0) {
	    std::swap(tets[i].vi[2], tets[i].vi[3]);
	}
    }

    cerr<<"min value: "<<min_value<<endl;
    cerr<<"max value: "<<max_value<<endl;

    fclose(f);
    return true;
}


void TetMesh::GetShell(TriangleMesh &tm, vector<real_type> &scalars) const {


    std::map< std::pair<int, std::pair<int,int> >, int> facekey3_to_index;	// facekey3 = sorted vertex indices
    vector< vector< std::pair<int,int> > > index_to_facekey2; // facekey2 = (tetindex, faceindex)

    for (unsigned t=0; t<tets.size(); t++) {
	for (int f=0; f<4; f++) {
	    int v1,v2,v3;

	    switch (f) {
	    case 0:
		v1=0; v2=1; v3=2;	break;
	    case 1:
		v1=0; v2=2; v3=3;	break;
	    case 2:
		v1=0; v2=3; v3=1;	break;
	    case 3:
		v1=1; v2=3; v3=2;	break;
	    }

	    v1 = tets[t].vi[v1];
	    v2 = tets[t].vi[v2];
	    v3 = tets[t].vi[v3];

	    // sort the vert indices
	    if (v2<v1)	std::swap(v1,v2);
	    if (v3<v2)	std::swap(v2,v3);
	    if (v2<v1)	std::swap(v1,v2);


	    std::pair<int, std::pair<int,int> > fkey3(v1, std::pair<int,int>(v2,v3));


	    std::map< std::pair<int, std::pair<int,int> >, int>::iterator it = facekey3_to_index.find(fkey3);


	    std::pair<int,int> fkey2(t,f);

	    if (it == facekey3_to_index.end()) {
		facekey3_to_index.insert(std::pair<std::pair<int, std::pair<int,int> >, int>(fkey3, index_to_facekey2.size()));

		index_to_facekey2.resize(index_to_facekey2.size()+1);
		index_to_facekey2.back().push_back(fkey2);
	    } else {
		index_to_facekey2[it->second].push_back(fkey2);
	    }

	}

    }


    // add all the faces that only appeared once
    IndexedTriangleSet its;
    std::map<int,int> tet_to_tri_verts;
    scalars.clear();

    for (unsigned i=0; i<index_to_facekey2.size(); i++) {

	if (index_to_facekey2[i].size() > 2) {
	    cerr<<index_to_facekey2[i].size()<<" faces with same index?!?"<<endl;
	}

	if (index_to_facekey2[i].size() == 1) {

	    int v1,v2,v3;

	    switch (index_to_facekey2[i][0].second) {
	    case 0:
		v1=0; v2=1; v3=2;	break;
	    case 1:
		v1=0; v2=2; v3=3;	break;
	    case 2:
		v1=0; v2=3; v3=1;	break;
	    case 3:
		v1=1; v2=3; v3=2;	break;
	    }

	    int v[3] = { tets[index_to_facekey2[i][0].first].vi[v1], 
			 tets[index_to_facekey2[i][0].first].vi[v2],
			 tets[index_to_facekey2[i][0].first].vi[v3] };

	    int itsvi[3];
	    for (int j=0; j<3; j++) {
		std::map<int,int>::iterator it = tet_to_tri_verts.find(v[j]);
		if (it == tet_to_tri_verts.end()) {
		    itsvi[j] = its.num_vertices();
		    tet_to_tri_verts.insert(std::pair<int,int>(v[j],itsvi[j]));
		    its.add_vertex(verts[v[j]].point);
		    scalars.push_back(verts[v[j]].scalar);
		} else {
		    itsvi[j] = it->second;
		}
	    }

	    its.insert_triangle((unsigned)itsvi[0], (unsigned)itsvi[1], (unsigned)itsvi[2], false, 0.0);
	}
    }

    ITS2TM(its, tm);
}


void TetMesh::SmoothScalars(int passes) {

    vector< std::set<int> > onerings(verts.size());	// 1-rings of all the vertices, plus itself
  
    vector< vector<int> > vtets(verts.size());
    for (unsigned t=0; t<tets.size(); t++) {
	for (int vi=0; vi<4; vi++) {
	    vtets[tets[t].vi[vi]].push_back(t);
	}
    }


    for (unsigned v=0; v<verts.size(); v++) {
	for (unsigned ti=0; ti<vtets[v].size(); ti++) {
	    int t = vtets[v][ti];
	    for (int vi=0; vi<4; vi++) {
		if (onerings[v].find(tets[t].vi[vi]) == onerings[v].end()) {
		    onerings[v].insert(tets[t].vi[vi]);
		}
	    }
	}
    }
  
  
  
    // smooth the scalars
    for (int pass=0; pass<passes; pass++) {
	vector<real_type> sscalars(verts.size(), 0);
	vector<int> valences(sscalars.size(), 0);	// not really valence
    
	for (unsigned t=0; t<tets.size(); t++) {
	    for (int vi1=0; vi1<4; vi1++) {
		for (int vi2=0; vi2<4; vi2++) {
		    int v1 = tets[t].vi[vi1];
		    int v2 = tets[t].vi[vi2];
	  
		    sscalars[v1] += verts[v2].scalar;
		    valences[v1]++;
		}
	    }
	}
	for (unsigned v=0; v<verts.size(); v++) {
	    verts[v].scalar = sscalars[v] / (real_type)valences[v];
	}
    }
}


real_type TetMesh::radius_ratio(int t) const {

    const Point3 *points[4] = { &verts[tets[t].vi[0]].point,
				&verts[tets[t].vi[1]].point,
				&verts[tets[t].vi[2]].point,
				&verts[tets[t].vi[3]].point };
    
    real_type a = Point3::distance(*points[0], *points[3]);
    real_type b = Point3::distance(*points[0], *points[2]);
    real_type c = Point3::distance(*points[2], *points[3]);
    real_type d = Point3::distance(*points[1], *points[2]);
    real_type e = Point3::distance(*points[1], *points[3]);
    real_type f = Point3::distance(*points[0], *points[1]);

    real_type V = (*points[2]-*points[0]).cross(*points[1]-*points[0]).dot(*points[3]-*points[0]) / (2*3);
    real_type S = ((*points[1]-*points[0]).cross(*points[2]-*points[0]).length() +
		   (*points[1]-*points[0]).cross(*points[3]-*points[0]).length() +
		   (*points[2]-*points[0]).cross(*points[3]-*points[0]).length() +
		   (*points[2]-*points[1]).cross(*points[3]-*points[1]).length()) * 0.5;
    V = fabs(V);
    S = fabs(S);

    real_type incircle_rad = 3 * V / S;
    real_type s = (a*d + b*e + c*f)/2;
    real_type circum_rad_inv = 6 * V / sqrt(s*(s-a*d)*(s-b*e)*(s-c*f));

    return (incircle_rad * circum_rad_inv);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TetMeshProjector::TetMeshProjector(const TetMesh &m, real_type iso) : mesh(m), isovalue(iso) {
    // setup the kdtree
    kdtree = new kdtree_type();
    for (unsigned t=0; t<mesh.tets.size(); t++) {
        kdtree->Insert(*this, t);
    }
}


Box3 TetMeshProjector::bounding_box(int tet) const {

    Box3 box = Box3::bounding_box(mesh.verts[mesh.tets[tet].vi[0]].point, 
				  mesh.verts[mesh.tets[tet].vi[1]].point, 
				  mesh.verts[mesh.tets[tet].vi[2]].point);
    box.update(mesh.verts[mesh.tets[tet].vi[3]].point);
    return box;
}


real_type TetMeshProjector::distance(int tet, const Point3 &from) const {

    if (PointInTet(from, tet))
	return 0;


    Triangle3 tris[4] = { Triangle3(mesh.verts[mesh.tets[tet].vi[0]].point,
				    mesh.verts[mesh.tets[tet].vi[1]].point,
				    mesh.verts[mesh.tets[tet].vi[2]].point),
			  Triangle3(mesh.verts[mesh.tets[tet].vi[0]].point,
				    mesh.verts[mesh.tets[tet].vi[2]].point,
				    mesh.verts[mesh.tets[tet].vi[3]].point),
			  Triangle3(mesh.verts[mesh.tets[tet].vi[0]].point,
				    mesh.verts[mesh.tets[tet].vi[3]].point,
				    mesh.verts[mesh.tets[tet].vi[1]].point),
			  Triangle3(mesh.verts[mesh.tets[tet].vi[3]].point,
				    mesh.verts[mesh.tets[tet].vi[2]].point,
				    mesh.verts[mesh.tets[tet].vi[1]].point) };


    real_type dist = std::min(std::min(Point3::distance(from, tris[0].closest_point(from)),
				       Point3::distance(from, tris[1].closest_point(from))),
			      std::min(Point3::distance(from, tris[2].closest_point(from)),
				       Point3::distance(from, tris[3].closest_point(from))));
    dist -= 1e-5;
    return dist;
}


bool TetMeshProjector::PointInTet(const Point3 &p, int t) const {

    const Point3 &p0 = mesh.verts[mesh.tets[t].vi[0]].point;
    const Point3 &p1 = mesh.verts[mesh.tets[t].vi[1]].point;
    const Point3 &p2 = mesh.verts[mesh.tets[t].vi[2]].point;
    const Point3 &p3 = mesh.verts[mesh.tets[t].vi[3]].point;

    if ((p1-p0).cross(p3-p0).dot(p-p0) < 0) return false;
    if ((p2-p0).cross(p1-p0).dot(p-p0) < 0) return false;
    if ((p3-p0).cross(p2-p0).dot(p-p0) < 0) return false;
    if ((p2-p1).cross(p3-p1).dot(p-p1) < 0) return false;
    return true;
}


void TetMeshProjector::TetBaryCoords(int t, const Point3 &x, real_type b[4]) const {

    Point3 p[4] = {
	mesh.verts[mesh.tets[t].vi[0]].point,
	mesh.verts[mesh.tets[t].vi[1]].point,
	mesh.verts[mesh.tets[t].vi[2]].point,
	mesh.verts[mesh.tets[t].vi[3]].point,
    };

    real_type area = (p[2]-p[0]).cross(p[1]-p[0]).dot(p[3]-p[0]);
    if (area == 0) {
	b[0] = 0.25;
	b[1] = 0.25;
	b[2] = 0.25;
	b[3] = 0.25;
    } else {
	real_type areai = 1 / area;
	b[0] = gsub2(p[2],x).cross(gsub2(p[1],x)).dot(gsub2(p[3],x)) * areai;
	b[1] = gsub2(p[0],x).cross(gsub2(p[2],x)).dot(gsub2(p[3],x)) * areai;
	b[2] = gsub2(p[0],x).cross(gsub2(p[3],x)).dot(gsub2(p[1],x)) * areai;
	b[3] = gsub2(p[0],x).cross(gsub2(p[1],x)).dot(gsub2(p[2],x)) * areai;
    }
}

real_type TetMeshProjector::TetBaryInterp(int t, const Point3 &x, real_type f0, real_type f1, real_type f2, real_type f3) const {
    real_type b[4];
    TetBaryCoords(t, x, b);
    return (b[0]*f0 + b[1]*f1 + b[2]*f2 + b[3]*f3);
}




int TetMeshProjector::ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const {

#ifdef TET_CIRCULAR_PROJECTION


    vector<int> nbrs;
    StartRingProjection(fp, nbrs);

    Vector3 axis = base2.position - base1.position;
    axis.normalize();

    Point3 center = base1.position + (axis.dot(fp - base1.position)) * axis;
    real_type rad = Point3::distance(center, fp);


    Vector3 udir = fp - center;
    udir.normalize();
    Vector3 vdir = axis.cross(udir);
    vdir.normalize();

    udir *= rad;
    vdir *= rad;

    RingEval re(*this, center, udir, vdir, nbrs);


    // look for a decent bracket
    real_type theta_min, theta_max, f_theta_min, f_theta_max;
    theta_min = theta_max = 0;
    f_theta_min = f_theta_max = re(0);

    if (f_theta_max < 0) { 
		
	// turn towards positive
	do {

	    theta_min = theta_max;
	    f_theta_min = f_theta_max;

	    theta_max += theta_step;
	    f_theta_max = re(theta_max);

	    if (re.hit_boundary) {
		cerr<<"couldn't eval ring position"<<endl;
		return PROJECT_FAILURE;// should be boundary, but since we start at the boundary we dont expect that to happen
	    }

	    if (theta_max>M_PI_2) {
		cerr<<"couldn't bracket"<<endl;
		return PROJECT_FAILURE;
	    }

	} while (f_theta_max<0);

    } else {

	// turn towards negative
	do {

	    theta_max = theta_min;
	    f_theta_max = f_theta_min;

	    theta_min -= theta_step;
	    f_theta_min = re(theta_min);

	    if (re.hit_boundary) {
		cerr<<"couldn't eval ring position"<<endl;
		return PROJECT_FAILURE;
	    }

	    if (theta_min<-M_PI_2) {
		cerr<<"couldn't bracket"<<endl;
		return PROJECT_FAILURE;
	    }

	} while (f_theta_min>0);

    }

    real_type zero = NR::zbrent(re, theta_min, theta_max, zbrent_tol*Point3::distance(base1.position, base2.position));
    tp = center + (real_type)cos(zero)*udir + (real_type)sin(zero)*vdir;
    NormalAtPoint(tp, tn);
    return PROJECT_SUCCESS;


#else
    return ProjectPoint(fp, tp, tn);
#endif
}



int TetMeshProjector::GetPointTet(const Point3 &p) const {

    vector<int> possibles;
    kdtree->GetIntersectedBoxes(*this, p, possibles);

    for (unsigned i=0; i<possibles.size(); i++) {
	if (PointInTet(p, possibles[i]))
	    return possibles[i];
    }


    if (allow_outside) {
	int ret;
	kdtree_type::OrderedTraverse ot(*kdtree, p, *this);
	ot.next(ret);
	return ret;
    }

    return -1;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////








template <typename T>
    void TetMeshProjectorNielson::HermiteCoefs(const T &f0, const T &d0, const T &fs, const T &ds, const T &s, T coefs[4]) {
    // f0 = f(0), fs=f(s) ...

    T ss = s*s;
    T sss = s*ss;

    coefs[0] = 2*f0/sss;
    coefs[1] =-3*f0/ss;
    coefs[2] = 0;
    coefs[3] = 1*f0;

    coefs[0] +=-2*fs/sss;
    coefs[1] += 3*fs/ss;

    coefs[0] += 1*d0/ss;
    coefs[1] +=-2*d0/s;
    coefs[2] += 1*d0;

    coefs[0] += 1*ds/ss;
    coefs[1] +=-1*ds/s;

}

template <typename T, typename R>
    void TetMeshProjectorNielson::HermiteEval(const T &f0, const T &d0, const T &fs, const T &ds, const T &s, const R &x, R &f) {
    T coefs[4];
    HermiteCoefs(f0, d0, fs, ds, s, coefs);
    f = ((coefs[0]*x + coefs[1])*x + coefs[2])*x + coefs[3];
}


template <typename T>
    void TetMeshProjectorNielson::EvalOnFace(const Point2 p[3],
					     const Vector2 g[3],
					     const real_type f[3],
					     const gtb::tPoint2<T> &x,
					     T &fx) const {

    T gfv[3];// = {	(T)f[0], (T)f[1], (T)f[2] };
    gfv[0] = f[0];
    gfv[1] = f[1];
    gfv[2] = f[2];
    const gtb::tPoint2<T> gpv[3]  = { gp<T>(p[0]), gp<T>(p[1]), gp<T>(p[2]) };
    const gtb::tVector2<T> ggv[3] = { gv<T>(g[0]), gv<T>(g[1]), gv<T>(g[2]) };



    real_type areai = 1 / cross(p[1]-p[0], p[2]-p[0]);
    const T gb[3] = { cross(gsub2(p[(0+1)%3],x), gsub2(p[(0+2)%3],x)) * areai,
		      cross(gsub2(p[(1+1)%3],x), gsub2(p[(1+2)%3],x)) * areai,
		      cross(gsub2(p[(2+1)%3],x), gsub2(p[(2+2)%3],x)) * areai };

    fx = 0;
    T wsum; wsum=0;
    for (int i=0; i<3; i++) {

	Vector2 edir = p[(i+2)%3] - p[(i+1)%3];
	real_type elen = edir.length();
	edir /= elen;

	T alpha = gb[(i+1)%3] / (gb[(i+1)%3] + gb[(i+2)%3]);
	T s = (1-alpha) * elen;

	const real_type ef0 = f[(i+1)%3];
	const real_type ed0 = edir.dot(g[(i+1)%3]);
	const real_type efs = f[(i+2)%3];
	const real_type eds = edir.dot(g[(i+2)%3]);

	Gradient<T,1> fsi;
	Gradient<T,1> gs = s;
	gs.gradient[0] = 1;
	HermiteEval(ef0, ed0, efs, eds, elen, gs, fsi);

	gtb::tVector2<T> dfsi = gscale(fsi.gradient[0], edir);	// along edge
	Vector2 enorm = edir.normal();	// perp to edge
	dfsi += gscale(gdot2(enorm, (gscale(alpha,g[(i+1)%3]) + gscale(1-alpha,g[(i+2)%3]))), enorm);

	gtb::tPoint2<T> si = gadd(p[(i+1)%3], gscale(s,edir));
	gtb::tVector2<T> dir = gsub1(si, p[i]);
	T max = dir.length();
	dir /= max;

	T f0; f0 = f[i];
	T d0 = gdot1(dir, g[i]);
	T fs = fsi.value;
	T ds = dir.dot(dfsi);

	T sx;
	HermiteEval(f0, d0, fs, ds, max, dir.dot(x-gpv[i]), sx);


	T w = gb[(i+1)%3] * gb[(i+2)%3];	// square this?
	fx += w*sx;
	wsum += w;

    }

    fx /= wsum;
}




template <typename T>
    void TetMeshProjectorNielson::EvalInTet(const Point3 p[4],
					    const Vector3 g[4],
					    const real_type f[4],
					    const gtb::tPoint3<T> &x,
					    T &fx) const {

    // assumption on ordering of vertices!
    real_type areai = 1 / (p[2]-p[0]).cross(p[1]-p[0]).dot(p[3]-p[0]);
    T gb[4] = { gsub2(p[2],x).cross(gsub2(p[1],x)).dot(gsub2(p[3],x)) * areai,
		gsub2(p[0],x).cross(gsub2(p[2],x)).dot(gsub2(p[3],x)) * areai,
		gsub2(p[0],x).cross(gsub2(p[3],x)).dot(gsub2(p[1],x)) * areai,
		gsub2(p[0],x).cross(gsub2(p[1],x)).dot(gsub2(p[2],x)) * areai };

    T ming = std::min(std::min(gb[0],gb[1]), std::min(gb[2],gb[3]));

    if (1||ming < 1e-4) {
	//	if (false) {//ming < 1e-2) {

	T muck = (T)1e-4; //2*ming;
	if (muck<0)
	    muck = -muck;

	muck += (T)1e-3;

	T nsum = (T)0;
	for (int i=0; i<4; i++) {
	    gb[i] += muck;
	    nsum += gb[i];
	}
	for (int i=0; i<4; i++) {
	    gb[i] /= nsum;
	}
    }

    T nming = std::min(std::min(gb[0],gb[1]), std::min(gb[2],gb[3]));
    if (nming <= 0) {
	//	  cerr<<"min still less than zero!"<<endl;
    }


    fx = 0;
    T wsum=(T)0;
    for (int i=0; i<4; i++) {

	Vector3 udir = p[(i+2)%4] - p[(i+1)%4];
	Vector3 vdir = p[(i+3)%4] - p[(i+1)%4];
	Vector3 norm = udir.cross(vdir);
	vdir = norm.cross(udir);
	norm.normalize();
	udir.normalize();
	vdir.normalize();

	T bsumi = 1 / (gb[(i+1)%4] + gb[(i+2)%4] + gb[(i+3)%4]);
	gtb::tPoint3<T> si((T)0, (T)0, (T)0);
	gadd_scaled(si, p[(i+1)%4], gb[(i+1)%4] * bsumi);
	gadd_scaled(si, p[(i+2)%4], gb[(i+2)%4] * bsumi);
	gadd_scaled(si, p[(i+3)%4], gb[(i+3)%4] * bsumi);

	gtb::tPoint2<T> s2(gdot2(udir,gsub1(si,p[(i+1)%4])), gdot2(vdir,gsub1(si,p[(i+1)%4])));

	Point2 p2[3] = { Point2(udir.dot(p[(i+1)%4]-p[(i+1)%4]), vdir.dot(p[(i+1)%4]-p[(i+1)%4])),
			 Point2(udir.dot(p[(i+2)%4]-p[(i+1)%4]), vdir.dot(p[(i+2)%4]-p[(i+1)%4])),
			 Point2(udir.dot(p[(i+3)%4]-p[(i+1)%4]), vdir.dot(p[(i+3)%4]-p[(i+1)%4])) };
	Vector2 g2[3] = { Vector2(udir.dot(g[(i+1)%4]), vdir.dot(g[(i+1)%4])),
			  Vector2(udir.dot(g[(i+2)%4]), vdir.dot(g[(i+2)%4])),
			  Vector2(udir.dot(g[(i+3)%4]), vdir.dot(g[(i+3)%4])) };
	real_type f2[3] = { f[(i+1)%4], f[(i+2)%4], f[(i+3)%4] };


	Gradient<T,2> fsi;
	gtb::tPoint2< Gradient<T,2> > gsi(s2[0],s2[1]);
	gsi[0].gradient[0] = 1;
	gsi[1].gradient[1] = 1;
	EvalOnFace(p2, g2, f2, gsi, fsi);


	gtb::tVector3<T> dfsi = gscale(fsi.gradient[0],udir) + gscale(fsi.gradient[1],vdir);	// parallel to the face
	dfsi += gscale(gdot2(norm, (bsumi*(gscale(gb[(i+1)%4],g[(i+1)%4]) +
					   gscale(gb[(i+2)%4],g[(i+2)%4]) + 
					   gscale(gb[(i+3)%4],g[(i+3)%4])))), norm);		// perp to the face


	gtb::tVector3<T> dir = gsub1(si, p[i]);
	T max = dir.length();
	dir /= max;

	T f0 = (T)f[i];
	T d0 = gdot1(dir, g[i]);
	T fs = fsi.value;
	T ds = dir.dot(dfsi);

	T sx;
	HermiteEval(f0, d0, fs, ds, max, dir.dot(gsub1(x,p[i])), sx);


	T w = gb[(i+1)%4] * gb[(i+2)%4] * gb[(i+3)%4];		w *= w;
	fx += w*sx;
	wsum += w;

    }

    fx /= wsum;
}


template <typename T>
    void TetMeshProjectorNielson::EvalInTet(int t, const gtb::tPoint3<T> &x, T &fx) const {

    gtb::tPoint3<real_type> dp[4]	= { gp<real_type>(mesh.verts[mesh.tets[t].vi[0]].point),
					    gp<real_type>(mesh.verts[mesh.tets[t].vi[1]].point),
					    gp<real_type>(mesh.verts[mesh.tets[t].vi[2]].point),
					    gp<real_type>(mesh.verts[mesh.tets[t].vi[3]].point) };
    gtb::tVector3<real_type> dg[4]	= { gv<real_type>(gradients[mesh.tets[t].vi[0]]),
					    gv<real_type>(gradients[mesh.tets[t].vi[1]]),
					    gv<real_type>(gradients[mesh.tets[t].vi[2]]),
					    gv<real_type>(gradients[mesh.tets[t].vi[3]]) };
    real_type ds[4]					= { (real_type)(mesh.verts[mesh.tets[t].vi[0]].scalar),
							    (real_type)(mesh.verts[mesh.tets[t].vi[1]].scalar),
							    (real_type)(mesh.verts[mesh.tets[t].vi[2]].scalar),
							    (real_type)(mesh.verts[mesh.tets[t].vi[3]].scalar) };

    EvalInTet(dp, dg, ds, x, fx);
}


TetMeshProjectorNielson::TetMeshProjectorNielson(const TetMesh &m, real_type iso) : TetMeshProjector(m, iso) {

    gradients.resize(mesh.verts.size(), Vector3(0,0,0));

    vector< vector<int> > vtets(mesh.verts.size());
    for (unsigned t=0; t<mesh.tets.size(); t++) {
	for (int vi=0; vi<4; vi++) {
	    vtets[mesh.tets[t].vi[vi]].push_back(t);
	}
    }


    vector< std::set<int> > onerings(mesh.verts.size());	// 1-rings of all the vertices, plus itself

    for (unsigned v=0; v<mesh.verts.size(); v++) {
	for (unsigned ti=0; ti<vtets[v].size(); ti++) {
	    int t = vtets[v][ti];
	    for (int vi=0; vi<4; vi++) {
		if (onerings[v].find(mesh.tets[t].vi[vi]) == onerings[v].end()) {
		    onerings[v].insert(mesh.tets[t].vi[vi]);
		}
	    }
	}
    }




    for (unsigned v=0; v<mesh.verts.size(); v++) {

	// put the 1-ring into amat/avec so it can be solved by the poly class
	const int numrings = 2;
	vector< std::set<int> > rings(numrings+1);
	rings[0].insert(v);

	for (int r=0; r<numrings; r++) {
	    rings[r+1] = rings[r];
	    for (std::set<int>::iterator i=rings[r].begin(); i!=rings[r].end(); ++i) {
		rings[r+1].insert(onerings[*i].begin(), onerings[*i].end());
	    }
	}

	// sort the points and get the median distance
	vector< std::pair<real_type,int> > sorted;
	for (std::set<int>::iterator i=rings[numrings].begin(); i!=rings[numrings].end(); ++i) {
	    sorted.push_back(std::pair<real_type,int>(Point3::distance(mesh.verts[v].point,mesh.verts[*i].point), *i));
	}
	std::sort(sorted.begin(), sorted.end());
	real_type mediandist = sorted[sorted.size()/2].first;



	vector<real_type> xs, ys, zs, fxs, weights;
	for (unsigned r=0; r<rings.size(); r++) {
	    for (std::set<int>::iterator i=rings[r].begin(); i!=rings[r].end(); ++i) {
		xs.push_back(mesh.verts[*i].point[0]);
		ys.push_back(mesh.verts[*i].point[1]);
		zs.push_back(mesh.verts[*i].point[2]);
		fxs.push_back(mesh.verts[*i].scalar);
		weights.push_back(exp(-Point3::squared_distance(mesh.verts[v].point,mesh.verts[*i].point)/(mediandist*mediandist)));
		//				cerr<<"fixme"<<endl;
	    }
	}


	const int degree=2;
	gtb::Poly3<degree, real_type, gtb::Monomial> poly;
	poly.LeastSquares(xs, ys, zs, fxs, weights, lls_normal_equations , lls_solver);



	real_type x = mesh.verts[v].point[0];
	real_type y = mesh.verts[v].point[1];
	real_type z = mesh.verts[v].point[2];
	for (int ix=0; ix <= degree; ++ix) {
	    for (int iy = 0; iy <= degree - ix; ++iy) {
		for (int iz = 0; iz <= degree - (iy+ix); ++iz) {
		    if (ix>0)
			gradients[v][0] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    ix * poly.TermEval(ix-1, x) * poly.TermEval(iy, y) * poly.TermEval(iz, z);
		    if (iy>0)
			gradients[v][1] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iy * poly.TermEval(ix, x) * poly.TermEval(iy-1, y) * poly.TermEval(iz, z);
		    if (iz>0)
			gradients[v][2] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iz * poly.TermEval(ix, x) * poly.TermEval(iy, y) * poly.TermEval(iz-1, z);
		}
	    }
	}
    }


    // smooth the gradients
    for (int sp=0; sp<nielson_gradient_smoothpasses; sp++) {
	vector<Vector3> sgradients(gradients.size(), Vector3(0,0,0));
	vector<int> valences(gradients.size(), 0);	// not really valence

	for (unsigned t=0; t<mesh.tets.size(); t++) {
	    for (int vi1=0; vi1<4; vi1++) {
		for (int vi2=0; vi2<4; vi2++) {
		    int v1 = mesh.tets[t].vi[vi1];
		    int v2 = mesh.tets[t].vi[vi2];

		    sgradients[v1] += gradients[v2];
		    valences[v1]++;
		}
	    }
	}
	for (unsigned v=0; v<mesh.verts.size(); v++) {
	    gradients[v] = sgradients[v] / (real_type)valences[v];
	}
    }
}

TetMeshProjectorNielson::~TetMeshProjectorNielson() {
    delete kdtree;
}




int TetMeshProjectorNielson::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const {

    // newton step onto the surface
    tp = fp;
    int cur_tet = -1;


    int iter=0;
    while (1) {

	iter++;
	if (iter>50) {
	    cerr<<"newton iteration failed to converge"<<endl;
	    return PROJECT_FAILURE;
	}

	if (cur_tet<0 || !PointInTet(tp, cur_tet)) {
	    cur_tet = GetPointTet(tp);
	}

	if (cur_tet < 0)	// outside the volume
	    return PROJECT_BOUNDARY;


	Gradient<projection_type,3> eval;
	gtb::tPoint3< Gradient<projection_type,3> > gx = gp< Gradient<projection_type,3> >(tp);
	gx[0].gradient[0] = 1;
	gx[1].gradient[1] = 1;
	gx[2].gradient[2] = 1;
	EvalInTet(cur_tet, gx, eval);

	projection_type step = -(eval.value-isovalue) / (eval.gradient[0]*eval.gradient[0] + eval.gradient[1]*eval.gradient[1] + eval.gradient[2]*eval.gradient[2]);

	tp[0] += (real_type)(step*eval.gradient[0]);
	tp[1] += (real_type)(step*eval.gradient[1]);
	tp[2] += (real_type)(step*eval.gradient[2]);

	tn[0]=(real_type)eval.gradient[0];
	tn[1]=(real_type)eval.gradient[1];
	tn[2]=(real_type)eval.gradient[2];


	// check for stopping
	if (step<step_dist_tol && fabs(eval.value-isovalue)<iso_dist_tol) {
	    break;
	}

    }

    tn.normalize();

    return PROJECT_SUCCESS;
									
}


bool TetMeshProjectorNielson::EvalAtPoint(const Point3 &p, real_type &f, vector<int> *nbrs) const {

    int t = GetPointTet(p);
    if (t<0) return false;

    // check for evaluating exactly at a mesh vertex
    for (int i=0; i<4; i++) {
	const Point3 &mp = mesh.verts[mesh.tets[t].vi[i]].point;
	if (p[0]==mp[0] && p[1]==mp[1] && p[2]==mp[2]) {
	    f = mesh.verts[mesh.tets[t].vi[i]].scalar;
	    return true;
	}
    }

    EvalInTet(t, p, f);
    return true;
}


bool TetMeshProjectorNielson::NormalAtPoint(const Point3 &p, Vector3 &n) const {

    int t = GetPointTet(p);
    if (t<0) return false;

    // check for evaluating exactly at a mesh vertex
    for (int i=0; i<4; i++) {
	const Point3 &mp = mesh.verts[mesh.tets[t].vi[i]].point;
	if (p[0]==mp[0] && p[1]==mp[1] && p[2]==mp[2]) {
	    n = gradients[mesh.tets[t].vi[i]];
	    n.normalize();
	    return true;
	}
    }

    gtb::tPoint3< Gradient<normal_type,3> > pg = gp< Gradient<normal_type,3> >(p);
    pg[0].gradient[0] = 1;
    pg[1].gradient[1] = 1;
    pg[2].gradient[2] = 1;

    Gradient<normal_type,3> fx;
    EvalInTet(t, pg, fx);

    n = Vector3(fx.gradient[0], fx.gradient[1], fx.gradient[2]);

    if (n[0] != n[0] ||
	n[1] != n[1] ||
	n[2] != n[2]) {
	n = gradients[mesh.tets[t].vi[0]];
    }

    n.normalize();
    return true;
}



bool TetMeshProjectorNielson::CurvatureAtPoint(const Point3 &p, real_type &k1, real_type &k2) const {

    int t = GetPointTet(p);
    if (t<0) return false;


    // check for evaluating exactly at a mesh vertex
    for (int i=0; i<4; i++) {
	const Point3 &mp = mesh.verts[mesh.tets[t].vi[i]].point;
	if (p[0]==mp[0] && p[1]==mp[1] && p[2]==mp[2]) {
	    return false;
	}
    }


    gtb::tPoint3< Gradient< Gradient<curvature_type,3> ,3> > pg = gp< Gradient< Gradient<curvature_type,3> ,3> >(p);
    pg[0].gradient[0].value = 1;
    pg[1].gradient[1].value = 1;
    pg[2].gradient[2].value = 1;
    pg[0].value.gradient[0] = 1;
    pg[1].value.gradient[1] = 1;
    pg[2].value.gradient[2] = 1;


    Gradient< Gradient<curvature_type,3> ,3> fx;
    EvalInTet(t, pg, fx);


    typedef gtb::tmat3<curvature_comp_type> mat3;
    mat3 H;
    H[0][0] = fx.gradient[0].gradient[0];
    H[0][1] = fx.gradient[0].gradient[1];
    H[0][2] = fx.gradient[0].gradient[2];
    H[1][0] = fx.gradient[1].gradient[0];
    H[1][1] = fx.gradient[1].gradient[1];
    H[1][2] = fx.gradient[1].gradient[2];
    H[2][0] = fx.gradient[2].gradient[0];
    H[2][1] = fx.gradient[2].gradient[1];
    H[2][2] = fx.gradient[2].gradient[2];

    curvature_comp_type gradient[3];
    gradient[0] = fx.gradient[0].value;
    gradient[1] = fx.gradient[1].value;
    gradient[2] = fx.gradient[2].value;


    gtb::tVector3<curvature_comp_type> normal(gradient[0], gradient[1], gradient[2]);
    curvature_comp_type gradmag = normal.length();
    normal /= gradmag;

    mat3 nnt(normal[0]*normal, normal[1]*normal, normal[2]*normal);
    mat3 I(gtb::tVector3<curvature_comp_type>(1,0,0), gtb::tVector3<curvature_comp_type>(0,1,0), gtb::tVector3<curvature_comp_type>(0,0,1));
    mat3 P = I - nnt;

    mat3 G = (-1/gradmag) * P*H*P;


    curvature_comp_type T = G[0][0] + G[1][1] + G[2][2];	// trace
    curvature_comp_type F=0;	// frobenius norm
    for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	    F += G[i][j]*G[i][j];
	}
    }
    F = sqrt(F);

    curvature_comp_type st = sqrt(std::max(0.0, 2*F*F-T*T));
    k1 = (real_type)((T+st)/2);
    k2 = (real_type)((T-st)/2);

    return true;
}








////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TetMeshProjectorMLS::SetRadiusParallel(int nt, int id, CSObject &cs) {


    for (int v=id; v<kdGetPoint.NumPoints(); v+=nt) {

	Point3 p = kdGetPoint(v);
	kdtree_type::OrderedIncrementalTraverse kdOrderedTraverse(*kdtree, p);

	int last = -1;

	for (int i=0; i<desired_knn; i++) {
	    real_type sd;
	    int next = kdOrderedTraverse.GetNext(sd);
	    if (next < 0) break;
	    last = next;
	}

	if (last >= 0) {
	    if (v < (int)mesh.verts.size())
		vert_radius[v] = Point3::distance(p, kdGetPoint(last));
	    else
		kdGetPoint.extra_radius[v-mesh.verts.size()] = Point3::distance(p, kdGetPoint(last));
	}
	else
	    cerr<<"aoeu"<<endl;
    }



}


void testInsideOutside(int nt, int id, const TetMeshProjector &tmp, const TetMesh &mesh) {

    for (unsigned t=id; t<mesh.tets.size(); t+=nt) {

	for (int i=0; i<50; i++) {

      
	    real_type b0 = myran1f(id);
	    real_type b1 = myran1f(id);
	    real_type b2 = myran1f(id);
	    real_type b3 = myran1f(id);

	    real_type sum = b0+b1+b2+b3;
	    b0 /= sum;
	    b1 /= sum;
	    b2 /= sum;
	    b3 /= sum;


	    Point3 p(0,0,0);
	    p.add_scaled(mesh.verts[mesh.tets[t].vi[0]].point, b0);
	    p.add_scaled(mesh.verts[mesh.tets[t].vi[1]].point, b1);
	    p.add_scaled(mesh.verts[mesh.tets[t].vi[2]].point, b2);
	    p.add_scaled(mesh.verts[mesh.tets[t].vi[3]].point, b3);


      
	    bool outside;

	    if (i>0)    outside = tmp.GetPointTet(p) < 0;
	    else        outside = !tmp.PointInTet(p, t);

	    if (outside) {
		cerr<<"outside!"<<endl;

		static thlib::CSObject cs;
		cs.enter();
		DbgPoints::add(p, 1, 0, 0);
		DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0));
		DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0));
		DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0));
		DbgPoly::add(mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0),
			     mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0));
	
		redrawAndWait(' ');
		dbgClear();
		cs.leave();

	    }

      
      

	}

    }

}


TetMeshProjectorMLS::TetMeshProjectorMLS(const TetMesh &m, real_type iso) : TetMeshProjector(m, iso), kdGetPoint() {


    kdGetPoint.mlsproj = this;


    //	ParallelExecutor(idealNumThreads, testInsideOutside, *this, m);
    //	cerr<<"aoeu"<<endl;
    //	redrawAndWait(' ');



    // setup the kdtree
    Box3 bbox = bounding_box(0);
    for (unsigned t=1; t<mesh.tets.size(); t++) {
	bbox = Box3::make_union(bbox, bounding_box(t));
    }



    for (unsigned t=0; t<mesh.tets.size(); t++) {

	/*
	  dbgClear();
	  DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0));
	  DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0));
	  DbgPoly::add(mesh.verts[mesh.tets[t].vi[0]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0));
	  DbgPoly::add(mesh.verts[mesh.tets[t].vi[1]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[2]].point, Point3(0,1,0), 
	  mesh.verts[mesh.tets[t].vi[3]].point, Point3(0,1,0));
	  dbgClear();
	  DbgPoints::add(mesh.verts[mesh.tets[t].vi[0]].point, 0,1,0);
	  DbgPoints::add(mesh.verts[mesh.tets[t].vi[1]].point, 0,1,0);
	  DbgPoints::add(mesh.verts[mesh.tets[t].vi[2]].point, 0,1,0);
	  DbgPoints::add(mesh.verts[mesh.tets[t].vi[3]].point, 0,1,0);
	*/

	// sample the interiors
	for (int b0=1; b0<eval_sub; b0++) {
	    for (int b1=1; b1<b0; b1++) {
		for (int b2=1; b2<b1; b2++) {
		
		    const int *vi = m.tets[t].vi;
		
		    real_type rb0, rb1, rb2, rb3;
	      
		    rb0 = (b2) / (real_type)eval_sub;
		    rb1 = (b1) / (real_type)eval_sub - rb0;
		    rb2 = (b0) / (real_type)eval_sub - rb0 - rb1;
		    rb3 = 1 - rb0 - rb1 - rb2;
		

		    Point3 p(0,0,0);
		    p.add_scaled(mesh.verts[vi[0]].point, rb0);
		    p.add_scaled(mesh.verts[vi[1]].point, rb1);
		    p.add_scaled(mesh.verts[vi[2]].point, rb2);
		    p.add_scaled(mesh.verts[vi[3]].point, rb3);
	      
		    real_type s=0;
		    s += mesh.verts[vi[0]].scalar * rb0;
		    s += mesh.verts[vi[1]].scalar * rb1;
		    s += mesh.verts[vi[2]].scalar * rb2;
		    s += mesh.verts[vi[3]].scalar * rb3;

		    kdGetPoint.extra_pts.push_back(p);
		    kdGetPoint.extra_scalars.push_back(s);

		    //		DbgPoints::add(p, 1,0,0);

		}
	    }
	}


	// sample the faces
	for (int f=0; f<4; f++) {

	    int vi[3];
	    switch (f) {
	    case 0: vi[0]=0; vi[1]=1; vi[2]=2; break;
	    case 1: vi[0]=0; vi[1]=2; vi[2]=3; break;
	    case 2: vi[0]=0; vi[1]=3; vi[2]=1; break;
	    case 3: vi[0]=1; vi[1]=3; vi[2]=2; break;
	    }



	    for (int i=0; i<3; i++) {
		vi[i] = mesh.tets[t].vi[vi[i]];
	    }
	    
	    int mini = (vi[0]<vi[1]) ? 0 : 1;
	    if (vi[2]<vi[mini]) mini=2;
	    
	    int maxi = (vi[0]>vi[1]) ? 0 : 1;
	    if (vi[2]>vi[maxi]) maxi=2;
	    

	    if (maxi == (mini+1)%3)
		continue;


	    for (int b0=1; b0<eval_sub; b0++) {
		for (int b1=1; b1<b0; b1++) {
		  
		    real_type rb0, rb1, rb2;
		    rb0 = (b1) / (real_type)eval_sub;
		    rb1 = (b0) / (real_type)eval_sub - rb0;
		    rb2 = 1 - rb0 - rb1;
		

		    Point3 p(0,0,0);
		    p.add_scaled(mesh.verts[vi[0]].point, rb0);
		    p.add_scaled(mesh.verts[vi[1]].point, rb1);
		    p.add_scaled(mesh.verts[vi[2]].point, rb2);

	      
		    real_type s=0;
		    s += mesh.verts[vi[0]].scalar * rb0;
		    s += mesh.verts[vi[1]].scalar * rb1;
		    s += mesh.verts[vi[2]].scalar * rb2;

		    kdGetPoint.extra_pts.push_back(p);
		    kdGetPoint.extra_scalars.push_back(s);

		    //		DbgPoints::add(p, 0,0,1);

		}
	    }
	}

	//	  redrawAndWait(' ');
    }
	  

    // sample the edges
    // make a huge set of all the edges
    std::set< std::pair<int,int> > all_edges;
    for (unsigned t=0; t<m.tets.size(); t++) {
	for (int e=0; e<6; e++) {
	    int vi[2];
	    switch (e) {
	    case 0: vi[0]=0; vi[1]=1; break;
	    case 1: vi[0]=0; vi[1]=2; break;
	    case 2: vi[0]=0; vi[1]=3; break;
	    case 3: vi[0]=1; vi[1]=2; break;
	    case 4: vi[0]=2; vi[1]=3; break;
	    case 5: vi[0]=3; vi[1]=1; break;
	    }

	    vi[0] = m.tets[t].vi[vi[0]];
	    vi[1] = m.tets[t].vi[vi[1]];

	    if (vi[0] > vi[1])
		std::swap(vi[0], vi[1]);

	    all_edges.insert(std::pair<int,int>(vi[0],vi[1]));
	}
    }

    for (std::set< std::pair<int,int> >::iterator i=all_edges.begin(); i!=all_edges.end(); ++i) {

	int vi[2] = { i->first, i->second };

	for (int b0=1; b0<eval_sub; b0++) {
	    real_type rb0, rb1;
	    rb0 = b0 / (real_type)eval_sub;
	    rb1 = 1 - rb0;

	    Point3 p(0,0,0);
	    p.add_scaled(mesh.verts[vi[0]].point, rb0);
	    p.add_scaled(mesh.verts[vi[1]].point, rb1);

	      
	    real_type s=0;
	    s += mesh.verts[vi[0]].scalar * rb0;
	    s += mesh.verts[vi[1]].scalar * rb1;

	    kdGetPoint.extra_pts.push_back(p);
	    kdGetPoint.extra_scalars.push_back(s);
	    
	}
    }


    kdtree = new kdtree_type(10, bbox, kdGetPoint);
    for (unsigned i=0; i<(mesh.verts.size()+kdGetPoint.extra_pts.size()); i++)
        kdtree->Insert(i);
    kdtree->MakeTree();


    // get the average edge lengths to the 1-ring nbrs
    vector< vector<int> > vtets(mesh.verts.size());
    for (unsigned t=0; t<mesh.tets.size(); t++) {
	for (int vi=0; vi<4; vi++) {
	    vtets[mesh.tets[t].vi[vi]].push_back(t);
	}
    }


    vert_radius.resize(mesh.verts.size());
    kdGetPoint.extra_radius.resize(kdGetPoint.extra_pts.size());

    cerr<<"setting vertex radii"<<endl;
    ParallelExecutor(idealNumThreads, makeClassFunctor(this, &TetMeshProjectorMLS::SetRadiusParallel));
    cerr<<"ok"<<endl;
}


TetMeshProjectorMLS::~TetMeshProjectorMLS() {
    delete kdtree;
}


template <typename T>
    bool TetMeshProjectorMLS::EvalAtPoint(const Point3 &p, T *f, gtb::tVector3<T> *gradient, gtb::tmat3<T> *hessian, const vector<int> *nbrs) const {


    int cur_tet = GetPointTet(p);
    if (cur_tet < 0) {
	cerr<<"point not in domain"<<endl;
	return false;
    }
    real_type rad = TetBaryInterp(cur_tet, p, vert_radius[mesh.tets[cur_tet].vi[0]], 
				  vert_radius[mesh.tets[cur_tet].vi[1]], 
				  vert_radius[mesh.tets[cur_tet].vi[2]], 
				  vert_radius[mesh.tets[cur_tet].vi[3]]);

    vector<int> lnbrs;
    if (!nbrs) {
	kdtree->Extract(p, max_query_knn, cache_overestimate*rad, std::back_inserter(lnbrs));
	nbrs = &lnbrs;
    }



    vector<T> xs, ys, zs, fxs, weights;
    for (unsigned i=0; i<nbrs->size(); i++) {
	Vector3 diff = kdGetPoint((*nbrs)[i]) - p;
	xs.push_back(diff[0]);
	ys.push_back(diff[1]);
	zs.push_back(diff[2]);
	fxs.push_back(kdGetPoint.scalar((*nbrs)[i]));
	weights.push_back(exp(-diff.squared_length()/(2*rad*rad*weight_scale*weight_scale)));
    }

    const int degree=2;
    gtb::Poly3<degree, T, gtb::Monomial> poly;

    if (gtb::Poly3<degree, T, gtb::Monomial>::CoefCount > nbrs->size()) {
	cerr<<"not enough points to solve polynomial!"<<endl;
	DbgPoints::add(p, 1, 0, 0);
	redrawAndWait(' ');
	return false;
    }

    //	cerr<<nbrs->size()<<endl;

    poly.LeastSquares(xs, ys, zs, fxs, weights, lls_normal_equations, lls_solver);



    if (f) {
	*f = poly.Eval(0, 0, 0);
    }


    if (gradient) {

	(*gradient)[0] = (*gradient)[1] = (*gradient)[2] = 0;

	for (int ix=0; ix <= degree; ++ix) {
	    for (int iy = 0; iy <= degree - ix; ++iy) {
		for (int iz = 0; iz <= degree - (iy+ix); ++iz) {
		    if (ix>=1)
			(*gradient)[0] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    ix * poly.TermEval(ix-1, 0) * poly.TermEval(iy, 0) * poly.TermEval(iz, 0);
		    if (iy>=1)
			(*gradient)[1] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iy * poly.TermEval(ix, 0) * poly.TermEval(iy-1, 0) * poly.TermEval(iz, 0);
		    if (iz>=1)
			(*gradient)[2] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iz * poly.TermEval(ix, 0) * poly.TermEval(iy, 0) * poly.TermEval(iz-1, 0);
		}
	    }
	}
    }


    if (hessian) {
	(*hessian)[0][0] = (*hessian)[0][1] = (*hessian)[0][2] = 0;
	(*hessian)[1][0] = (*hessian)[1][1] = (*hessian)[1][2] = 0;
	(*hessian)[2][0] = (*hessian)[2][1] = (*hessian)[2][2] = 0;

	for (int ix=0; ix <= degree; ++ix) {
	    for (int iy = 0; iy <= degree - ix; ++iy) {
		for (int iz = 0; iz <= degree - (iy+ix); ++iz) {
		    if (ix>=2)
			(*hessian)[0][0] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    ix * (ix-1) * poly.TermEval(ix-2, 0) * poly.TermEval(iy, 0) * poly.TermEval(iz, 0);
		    if (ix>=1 && iy>=1)
			(*hessian)[0][1] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    ix * iy * poly.TermEval(ix-1, 0) * poly.TermEval(iy-1, 0) * poly.TermEval(iz, 0);
		    if (ix>=1 && iz>=1)
			(*hessian)[0][2] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    ix * iz * poly.TermEval(ix-1, 0) * poly.TermEval(iy, 0) * poly.TermEval(iz-1, 0);
		    if (iy>=2)
			(*hessian)[1][1] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iy * (iy-1) * poly.TermEval(ix, 0) * poly.TermEval(iy-2, 0) * poly.TermEval(iz, 0);
		    if (iy>=1 && iz>=1)
			(*hessian)[1][2] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iy * iz * poly.TermEval(ix, 0) * poly.TermEval(iy-1, 0) * poly.TermEval(iz-1, 0);
		    if (iz>=2)
			(*hessian)[2][2] += poly._coefs[poly.imap(ix, iy, iz)] * 
			    iz * (iz-1) * poly.TermEval(ix, 0) * poly.TermEval(iy, 0) * poly.TermEval(iz-2, 0);
		}
	    }
	}

	(*hessian)[1][0] = (*hessian)[0][1];
	(*hessian)[2][0] = (*hessian)[0][2];
	(*hessian)[2][1] = (*hessian)[1][2];
    }


    return true;
}



int TetMeshProjectorMLS::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const {


    int cur_tet = -1;

    // cache the nearest neighbors
    real_type rad = kdGetPoint.radius(kdtree->FindMin(fp));
	


    vector<int> nbrs;
    kdtree->Extract(fp, max_query_knn, cache_overestimate*rad, std::back_inserter(nbrs));



    // newton step onto the surface
    tp = fp;

    int iter=0;
    while (1) {

	iter++;
	if (iter>50) {
	    cerr<<"newton iteration failed to converge"<<endl;
	    return PROJECT_FAILURE;
	}


	if (!allow_outside) {
	    if (cur_tet<0 || !PointInTet(tp, cur_tet)) {
		cur_tet = GetPointTet(tp);
	    }
	    if (cur_tet < 0)	// outside the volume
		return PROJECT_BOUNDARY;
	}


	projection_type f;
	gtb::tVector3<projection_type>  gradient;
	if (!EvalAtPoint<projection_type>(tp, &f, &gradient, NULL, &nbrs))
	    return PROJECT_FAILURE;

	projection_type step = -(f-isovalue) / gradient.squared_length();

	tp[0] += (real_type)(step*gradient[0]);
	tp[1] += (real_type)(step*gradient[1]);
	tp[2] += (real_type)(step*gradient[2]);

	tn[0]=(real_type)gradient[0];
	tn[1]=(real_type)gradient[1];
	tn[2]=(real_type)gradient[2];


	// check for stopping
	if (step<step_dist_tol && fabs(f-isovalue)<iso_dist_tol) {
	    break;
	}
    }

    tn.normalize();
    return PROJECT_SUCCESS;
}


bool TetMeshProjectorMLS::EvalAtPoint(const Point3 &p, real_type &f, vector<int> *nbrs) const {

    eval_type e;
    if (EvalAtPoint<eval_type>(p, &e, NULL, NULL, nbrs)) {
	f = (real_type)e;
	return true;
    }
    return false;
}


bool TetMeshProjectorMLS::NormalAtPoint(const Point3 &p, Vector3 &n) const {

    gtb::tVector3<normal_type> nn;
    bool ret = EvalAtPoint<normal_type>(p, NULL, &nn, NULL);
    nn.normalize();
    n = Vector3(nn[0], nn[1], nn[2]);
    return ret;
}


bool TetMeshProjectorMLS::CurvatureAtPoint(const Point3 &p, real_type &k1, real_type &k2) const {

    gtb::tVector3<curvature_type> gradiente;
    gtb::tmat3<curvature_type> He;

    if (!EvalAtPoint<curvature_type>(p, NULL, &gradiente, &He))
	return false;

    gtb::tVector3<curvature_comp_type> gradient(gradiente[0], gradiente[1], gradiente[2]);
    gtb::tmat3<curvature_comp_type> H;
    for (int r=0; r<3; r++) {
	for (int c=0; c<3; c++) {
	    H[r][c] = He[r][c];
	}
    }


    gtb::tVector3<curvature_comp_type> normal = gradient;
    curvature_comp_type gradmag = normal.length();
    normal /= gradmag;

    gtb::tmat3<curvature_comp_type> nnt(normal[0]*normal, normal[1]*normal, normal[2]*normal);
    gtb::tmat3<curvature_comp_type> I(gtb::tVector3<curvature_comp_type>(1,0,0), gtb::tVector3<curvature_comp_type>(0,1,0), gtb::tVector3<curvature_comp_type>(0,0,1));
    gtb::tmat3<curvature_comp_type> P = I - nnt;

    gtb::tmat3<curvature_comp_type> G = (-1/gradmag) * P*H*P;


    curvature_comp_type T = G[0][0] + G[1][1] + G[2][2];	// trace
    curvature_comp_type F=0;	// frobenius norm
    for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	    F += G[i][j]*G[i][j];
	}
    }
    F = sqrt(F);

    curvature_comp_type st = sqrt(std::max((curvature_comp_type)0, 2*F*F-T*T));
    k1 = (real_type)((T+st)/2);
    k2 = (real_type)((T-st)/2);

    return true;
}


void TetMeshProjectorMLS::StartRingProjection(const Point3 &p, vector<int> &nbrs) const {

    nbrs.clear();

    // cache the nearest neighbors
    real_type rad = kdGetPoint.radius(kdtree->FindMin(p));
    kdtree->Extract(p, max_query_knn, cache_overestimate*rad, std::back_inserter(nbrs));
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MeshPointScalarEvalParallel(int nt, int id, CSObject &cs, const TetMeshProjector &tmp, const TetMesh &mesh, vector<real_type> &scalars) {
    for (unsigned v=id; v<mesh.verts.size(); v+=nt) {
	if (tmp.GetProjectorType() == TetMeshProjector::TET_MESH_PROJECTOR_NIELSON || 
	    !tmp.EvalAtPoint(mesh.verts[v].point, scalars[v], NULL))

	    scalars[v] = mesh.verts[v].scalar;
    }
}


void TetMeshGuidanceField::BuildGuidanceParallel(int nt, int id, CSObject &cs, real_type isovalue, const TetMesh &mesh, const vector<real_type> &scalars) { 


    vector<Point3> mypoints;
    vector<Vector3> mynorms;
    vector<real_type> mycurvatures;


    for (unsigned t=id; t<mesh.tets.size(); t+=nt) {

	real_type max = -1e34;
	real_type min = 1e34;

	for (int vi=0; vi<4; vi++) {
	    real_type f = scalars[mesh.tets[t].vi[vi]];
	    max = std::max(max, f);
	    min = std::min(min, f);
	}

	if (min<=isovalue && max>=isovalue) {

	    for (int b0=0; b0<curvature_sub; b0++) {
		for (int b1=0; b1<=b0; b1++) {
		    for (int b2=0; b2<=b1; b2++) {

			real_type rb0, rb1, rb2, rb3;

			rb0 = (b2+myran1f(id)) / curvature_sub;
			rb1 = (b1+myran1f(id)) / curvature_sub - rb0;
			rb2 = (b0+myran1f(id)) / curvature_sub - rb0 - rb1;
			rb3 = 1 - rb0 - rb1 - rb2;

			Point3 x(0,0,0);
			x.add_scaled(mesh.verts[mesh.tets[t].vi[0]].point, rb0);
			x.add_scaled(mesh.verts[mesh.tets[t].vi[1]].point, rb1);
			x.add_scaled(mesh.verts[mesh.tets[t].vi[2]].point, rb2);
			x.add_scaled(mesh.verts[mesh.tets[t].vi[3]].point, rb3);


			Point3 tp;
			Vector3 tn;
			real_type k1,k2;

			if (projector.ProjectPoint(x, tp, tn)==PROJECT_SUCCESS &&
			    projector.CurvatureAtPoint(tp, k1, k2)) {

			    real_type kmax = k2;
			    if (fabs(k1)>fabs(k2))
				kmax = k1;

			    mypoints.push_back(tp);
			    mynorms.push_back(tn);
			    mycurvatures.push_back(kmax);

			    if ((t-id)%(nt*1000) == 0) {
				cs.enter();
							  
				for (unsigned i=0; i<mypoints.size(); i++) {
				    kdGetPoint.allpoints.push_back(mypoints[i]);
				    DbgOPoints::add(mypoints[i], mynorms[i], 0, 1, 0);
				}
				for (unsigned i=0; i<mycurvatures.size(); i++) {
				    ideal_length.push_back(MaxCurvatureToIdeal(mycurvatures[i]));
				}
							  
				cs.leave();

				mypoints.clear();
				mynorms.clear();
				mycurvatures.clear();

			    }

			}
		    }
		}
	    }
	}
    }

    cs.enter();

    for (unsigned i=0; i<mypoints.size(); i++) {
	kdGetPoint.allpoints.push_back(mypoints[i]);
	DbgOPoints::add(mypoints[i], mynorms[i], 0, 1, 0);
    }
    for (unsigned i=0; i<mycurvatures.size(); i++) {
	ideal_length.push_back(MaxCurvatureToIdeal(mycurvatures[i]));
    }

    cs.leave();

}


TetMeshGuidanceField::TetMeshGuidanceField(TetMeshProjector &proj, real_type rho, real_type min_step, real_type max_step, real_type reduction)
    : GuidanceField(rho, min_step, max_step, reduction), projector(proj), kdGetPoint(), kdOrderedTraverse(NULL) {


    real_type isovalue = projector.GetIsoValue();
    const TetMesh &mesh = projector.GetMesh();


    cerr<<"evaluating at all points"<<endl;
    vector<real_type> scalars(mesh.verts.size());
    allow_outside = true;
    ParallelExecutor(idealNumThreads, &MeshPointScalarEvalParallel, projector, mesh, scalars);
    allow_outside = false;

    cerr<<"finding guidance field points"<<endl;
    ParallelExecutor(idealNumThreads, makeClassFunctor(this, &TetMeshGuidanceField::BuildGuidanceParallel), isovalue, mesh, scalars);



    redrawAndWait(' ');

    // setup the kdtree
    Box3 bbox = projector.bounding_box(0);
    for (unsigned t=1; t<mesh.tets.size(); t++) {
	bbox = Box3::make_union(bbox, projector.bounding_box(t));
    }

    kdtree = new kdtree_type(10, bbox, kdGetPoint);
    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++)
	kdtree->Insert(i);
    kdtree->MakeTree();


    vector<int> marked(ideal_length.size()); // for removal
    Trimmer<GetPoint> trimmer(*this, kdGetPoint, ideal_length);
    trimmer.Trim(marked);

    vector<Point3> keepers;
    vector<real_type> keep_curv;
    for (unsigned i=0; i<marked.size(); i++) {
	if (!marked[i]) {
	    keepers.push_back(kdGetPoint.allpoints[i]);
	    keep_curv.push_back(ideal_length[i]);
	}
    }
    kdGetPoint.allpoints = keepers;
    ideal_length = keep_curv;


    delete kdtree;
    kdtree = new kdtree_type(10, bbox, kdGetPoint);
    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++)
	kdtree->Insert(i);
    kdtree->MakeTree();


    cerr<<"Guidance field with "<<kdGetPoint.allpoints.size()<<" points"<<endl;
}


TetMeshGuidanceField::~TetMeshGuidanceField() {
    if (kdOrderedTraverse) delete kdOrderedTraverse;
    delete kdtree;
}


void* TetMeshGuidanceField::OrderedPointTraverseStart(const Point3 &p) {
    return new kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
}

bool TetMeshGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = kdGetPoint(n);
    ideal = ideal_length[n];
    return true;
}

void TetMeshGuidanceField::OrderedPointTraverseEnd(void *ctx) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    delete kdOrderedTraverse;
}



const Point3& TetMeshGuidanceField::PointLocation(int i) const {
    return kdGetPoint(i);
}


void TetMeshGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad) {
    pts.resize(kdGetPoint.allpoints.size());
    norms.resize(kdGetPoint.allpoints.size());
    rad.resize(kdGetPoint.allpoints.size());

    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++) {

	pts[i] = kdGetPoint.allpoints[i];
	if (!projector.NormalAtPoint(pts[i], norms[i])) norms[i] = Vector3(0,0,0);
	rad[i] = ideal_length[i];

    }
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mt_add_tri(int v00, int v01, real_type d00, real_type d01, 
		int v10, int v11, real_type d10, real_type d11, 
		int v20, int v21, real_type d20, real_type d21, 
		TriangulatorController &tc, const TetMesh &mesh,
		int &nverts, int &ntris) {


    int edges[3][2] = { { v00, v01 },
			{ v10, v11 }, 
			{ v20, v21 } };	

    real_type dists[3][2] = { { d00, d01 },
			      { d10, d11 }, 
			      { d20, d21 } };	


    Point3 vp[3];

    for (int i=0; i<3; i++) {

	real_type alpha = dists[i][0] / (dists[i][0]-dists[i][1]);

	vp[i] = mesh.verts[edges[i][0]].point;
	vp[i] += alpha * (mesh.verts[edges[i][1]].point - mesh.verts[edges[i][0]].point);
    }

    int vi[3];
    vi[0] = nverts++;
    vi[1] = nverts++;
    vi[2] = nverts++;

    Vector3 normal = (vp[1]-vp[0]).cross(vp[2]-vp[0]);
    normal.normalize();

    tc.AddVertex(vi[0], vp[0], normal, false);
    tc.AddVertex(vi[1], vp[1], normal, false);
    tc.AddVertex(vi[2], vp[2], normal, false);

    tc.AddTriangle(ntris++, vi[0], vi[1], vi[2]);
}


void MarchingTets(const TetMesh &mesh, TriangulatorController &tc, real_type isovalue) {

    int nverts=0;
    int ntris=0;


    for (unsigned t=0; t<mesh.tets.size(); t++) {

	int verts[4] = { mesh.tets[t].vi[0], mesh.tets[t].vi[1], mesh.tets[t].vi[2], mesh.tets[t].vi[3] };
	real_type dists[4] = { mesh.verts[verts[0]].scalar-isovalue,
			       mesh.verts[verts[1]].scalar-isovalue,
			       mesh.verts[verts[2]].scalar-isovalue,
			       mesh.verts[verts[3]].scalar-isovalue };

	bool matched = false;


	for (int v1=0; v1<4; v1++) {

	    if (dists[v1]>0 && dists[(v1+1)%4]<0 && dists[(v1+2)%4]<0 && dists[(v1+3)%4]<0) {

		mt_add_tri(verts[v1], verts[(v1+1)%4], dists[v1], dists[(v1+1)%4], 
			   verts[v1], verts[(v1+2)%4], dists[v1], dists[(v1+2)%4], 
			   verts[v1], verts[(v1+3)%4], dists[v1], dists[(v1+3)%4], 
			   tc, mesh, nverts, ntris);

		matched = true; break;
	    }


	    if (dists[v1]<0 && dists[(v1+1)%4]>0 && dists[(v1+2)%4]>0 && dists[(v1+3)%4]>0) {

		mt_add_tri(verts[v1], verts[(v1+1)%4], dists[v1], dists[(v1+1)%4], 
			   verts[v1], verts[(v1+3)%4], dists[v1], dists[(v1+3)%4], 
			   verts[v1], verts[(v1+2)%4], dists[v1], dists[(v1+2)%4], 
			   tc, mesh, nverts, ntris);

		matched = true; break;
	    }


	    for (int iv2=0; iv2<3; iv2++) {

		int v2 = (v1+1+(iv2+0)%3)%4;
		int v3 = (v1+1+(iv2+1)%3)%4;
		int v4 = (v1+1+(iv2+2)%3)%4;



		if (dists[v1]>0 && dists[v2]>0 && dists[v3]<0 && dists[v4]<0) {

		    mt_add_tri(verts[v2], verts[v3], dists[v2], dists[v3], 
			       verts[v1], verts[v3], dists[v1], dists[v3], 
			       verts[v1], verts[v4], dists[v1], dists[v4], 
			       tc, mesh, nverts, ntris);

		    mt_add_tri(verts[v1], verts[v4], dists[v1], dists[v4], 
			       verts[v2], verts[v4], dists[v2], dists[v4], 
			       verts[v2], verts[v3], dists[v2], dists[v3], 
			       tc, mesh, nverts, ntris);

		    matched = true; break;
		}



		if (dists[v1]<0 && dists[v2]<0 && dists[v3]>0 && dists[v4]>0) {

		    mt_add_tri(verts[v2], verts[v3], dists[v2], dists[v3], 
			       verts[v1], verts[v4], dists[v1], dists[v4], 
			       verts[v1], verts[v3], dists[v1], dists[v3], 
			       tc, mesh, nverts, ntris);

		    mt_add_tri(verts[v1], verts[v4], dists[v1], dists[v4], 
			       verts[v2], verts[v3], dists[v2], dists[v3], 
			       verts[v2], verts[v4], dists[v2], dists[v4], 
			       tc, mesh, nverts, ntris);

		    matched = true; break;
		}



	    }

	    if (matched) break;
	}


	if (matched) continue;


    }

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class TetLineEval {

    public:
    TetLineEval(const TetMeshProjector &t, const Point3 &_p0, const Point3 &_p1, real_type iso) :
		tetmeshproj(t), p0(_p0), p1(_p1), iso_value(iso) {
	}

    real_type operator()(real_type t) {

		Point3 x = p0 + t*(p1-p0);
		real_type v;
		tetmeshproj.EvalAtPoint(x, v, NULL);
		return v-iso_value;
    }

        
    const TetMeshProjector &tetmeshproj;
    real_type iso_value;
    Point3 p0;
    Point3 p1;
};


void MarchingTris(const TriangleMesh &mesh, const TetMeshProjector &tm, const vector<real_type> &scalar, real_type isovalue,
		  vector< vector<Point3> > &paths) {


	bool rootfind = true;

    paths.clear();

    vector<bool> donetris(mesh.faces.size(), false);


    int start=-1;

    while (1) {

	while (1) {
	    ++start;
	    if (start > (int)mesh.faces.size()) break;
	    if (!donetris[start]) {

		if (scalar[mesh.faces[start].verts[0]]-isovalue > 0 &&
		    scalar[mesh.faces[start].verts[1]]-isovalue > 0 &&
		    scalar[mesh.faces[start].verts[2]]-isovalue > 0)
		    donetris[start]=true;
		else if (scalar[mesh.faces[start].verts[0]]-isovalue < 0 &&
			 scalar[mesh.faces[start].verts[1]]-isovalue < 0 &&
			 scalar[mesh.faces[start].verts[2]]-isovalue < 0)
		    donetris[start]=true;
		else
		    break;
	    }
	}
	if (start > (int)mesh.faces.size()) break;

	int next = start;
	paths.push_back(vector<Point3>());

	do {

	    donetris[next]=true;

	    real_type dists[3] = { scalar[mesh.faces[next].verts[0]]-isovalue,
				   scalar[mesh.faces[next].verts[1]]-isovalue,
				   scalar[mesh.faces[next].verts[2]]-isovalue };

	    bool found=false;

	    for (int i=0; i<3; i++) {

		if (dists[i]<0 && dists[(i+1)%3]>0 && dists[(i+2)%3]>0) {
		    real_type alpha;
		    alpha = -dists[i] / (dists[(i+2)%3] - dists[i]);

			TetLineEval le(tm,
						   mesh.verts[mesh.faces[next].verts[i]].point, 
						   mesh.verts[mesh.faces[next].verts[(i+2)%3]].point, 
						   isovalue);

			if (rootfind) {
				extern real_type step_tol;
				extern real_type value_tol;
				real_type step = step_tol;

				for (int j=0; j<5; j++) {
					alpha = NR::zbrent(le, (real_type)0, (real_type)1, step);
					step *= 0.05;
					if (fabs(le(alpha)) < value_tol)
						break;
				}
			}

		    paths.back().push_back(mesh.verts[mesh.faces[next].verts[i]].point + 
					   alpha * (mesh.verts[mesh.faces[next].verts[(i+2)%3]].point - mesh.verts[mesh.faces[next].verts[i]].point));

		    next = mesh.faces[next].nbrs[i];
		    found=true; break;
		}

		if (dists[i]>0 && dists[(i+1)%3]<0 && dists[(i+2)%3]<0) {
		    real_type alpha;
		    alpha = -dists[i] / (dists[(i+1)%3] - dists[i]);

			TetLineEval le(tm,
						   mesh.verts[mesh.faces[next].verts[i]].point, 
						   mesh.verts[mesh.faces[next].verts[(i+1)%3]].point, 
						   isovalue);

			if (rootfind) {
				extern real_type step_tol;
				extern real_type value_tol;
				real_type step = step_tol;

				for (int j=0; j<5; j++) {
					alpha = NR::zbrent(le, (real_type)0, (real_type)1, step);
					step *= 0.05;
					if (fabs(le(alpha)) < value_tol)
						break;
				}
			}

		    paths.back().push_back(mesh.verts[mesh.faces[next].verts[i]].point + 
					   alpha * (mesh.verts[mesh.faces[next].verts[(i+1)%3]].point - mesh.verts[mesh.faces[next].verts[i]].point));
		    next = mesh.faces[next].nbrs[(i+2)%3];
		    found=true; break;
		}
	    }

	    if (!found) {
		cerr<<"error - didn't find isosurface crossing!"<<endl;
	    }

	} while (next != start);


    }


}




bool ProjectToIntersection(const SurfaceProjector  &surf1, const SurfaceProjector &surf2,
			   Point3 &p, Vector3 &norm1, Vector3 &norm2, real_type tol) {


    real_type move;

    int iter=0;
    do {
	iter++;
	if (iter > 100) {
	    cerr<<"too many iterations"<<endl;
	    allow_outside = false;
	    return false;
	}

	Point3 p1, p2;

	if (surf1.ProjectPoint(p, p1, norm1) != PROJECT_SUCCESS) {
	    cerr<<"projection1 failed"<<endl;
	    return false;
	}

	if (surf2.ProjectPoint(p1, p2, norm2) != PROJECT_SUCCESS) {
	    cerr<<"projection2 failed"<<endl;
	    return false;
	}

	move = std::max(Point3::distance(p,p1), Point3::distance(p,p2));
	p = p2;

    } while (move > tol);


    return true;
}







void WalkIntersectionCurve(const SurfaceProjector &p1, const SurfaceProjector &p2, 
			   const vector< vector<Point3> > &startpoints,
			   real_type stepsize,
			   vector< vector<Point3> > &out_points,
			   vector< vector<Vector3> > &out_norms1,
			   vector< vector<Vector3> > &out_norms2) {

    out_points.resize(startpoints.size());
    out_norms1.resize(startpoints.size());
    out_norms2.resize(startpoints.size());

    for (unsigned s=0; s<startpoints.size(); s++) {

	dbgClear();

	vector<Point3> &op = out_points[s];
	vector<Vector3> &on1 = out_norms1[s];
	vector<Vector3> &on2 = out_norms2[s];
 

	op.resize(1);
	on1.resize(1);
	on2.resize(1);


	// find the initial point
	const vector<Point3> &sp = startpoints[s];
	bool found=false;
	for (unsigned i=0; i<sp.size(); i++) {
     
	    op[0] = sp[i];
	    if (ProjectToIntersection(p1, p2, op[0], on1[0], on2[0], 1e-4)) {
		found=true;
		break;
	    }
	}

	if (!found) {
	    cerr<<"start not found!!!!"<<endl;
	    return;
	}


	DbgPoints::add(op[0], 1,0,0);

	// walk around the intersection with small steps until we get back to where we started
	while (1) {

	    Vector3 dir = on1.back().cross(on2.back());
	    dir.normalize();

	    Point3 nextp;
	    Vector3 nextn1;
	    Vector3 nextn2;
	
	    for (int attempt=0; attempt<1; attempt++) {

		nextp = op.back() + stepsize * dir;

		if (!ProjectToIntersection(p1, p2, nextp, nextn1, nextn2, 1e-4)) {
		    cerr<<"error stepping around the intersection curve!! smaller step needed?"<<endl;
		}


		if (1 || Point3::distance(op.back(), nextp) > stepsize*0.1) {
		    //	    (op.size()==0 || (op.back()-op[op.size()-2]).normalized().dot(nextp-op.back()) > stepsize*0.1)) {
		    op.push_back(nextp);
		    on1.push_back(nextn1);
		    on2.push_back(nextn2);
		    break;
		}
	    }

	    DbgPoints::add(op.back(), 0,1,0);
	    vector<Point3> line(2);
	    line[0] = op.back();
	    line[1] = line[0] + stepsize/2 * dir;
	    DbgPLines::add(line, 0);

	    if (op.size()>5 && Point3::distance(nextp,op[0])<2*stepsize)
		break;

	}
    }
}






