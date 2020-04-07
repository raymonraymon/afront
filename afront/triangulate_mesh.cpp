
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


//#define CLOSEST_POINT_PROJECTION


///////////////////////////////////////////////////////////////////////////////
// Helper for finding curvature of space curves (boundaries)

class Quadratic {
    public:

    Quadratic() {
	coefs[0] = 0;
	coefs[1] = 0;
	coefs[2] = 0;
    }

    Quadratic(real_type A, real_type B, real_type C) {
	coefs[0] = A;
	coefs[1] = B;
	coefs[2] = C;
    }

    real_type eval(real_type x) {
	return (coefs[0]*x*x + coefs[1]*x + coefs[2]);
    }

    real_type eval1(real_type x) {
	return (2*coefs[0]*x + coefs[1]);
    }

    real_type eval2(real_type x) {
	return (2*coefs[0]);
    }


    bool fit(const vector<real_type> &xi, const vector<real_type> &yi, const vector<real_type> &weights) {

	const int N = (int)xi.size();
	const int NCF = 3;
	const int DEG = 2;

	double cc[NCF];
	gtb::AMatc<double> A(NCF, NCF);
	gtb::AVec<double> b(NCF);
	A.set(0);
	b.set(0);

	real_type max_w = weights[0];
	for (int i=1; i<N; i++) {
	    if (weights[i] > max_w)
		max_w = weights[i];
	}
	double invmaxw = 1.0 / max_w;

	/*--------- Prepare the matrix A and vector b ----------*/
	for (int k = 0; k < N; ++k)
	    {
		double w = weights[k] * invmaxw;

		//
		// Speedup?
		//   - if our function is a*f_a(x,y) + b*f_b(x,y) + c*f_c(x,y) +...
		//     compute f_a(x,y), f_b(x,y), ...
		//     i.e. computes the "coefficients" of a,b,c...
		//   - Polynoms (old) pre-compute X^i, y^i, i=0:deg
		//
		{
		    int cidx = 0;
		    for (int p = 0; p <= DEG; ++p)
			{
			    cc[cidx] = gtb::ipow((double)xi[k], p);
			    ++cidx;
			}
		}

		double *pr, *pc;
		for (int i = 0; i < NCF; ++i)
		    {
			pc = pr = A.flatten() + i + i*NCF;
			for (int j = i; j < NCF; ++j)
			    {
				*pr = *pc += cc[i]*cc[j] *w;
				pr += NCF;
				++pc;
			    }

			// xi[2] is the height of the function above the 
			// plane. (i.e. y_i)
			b[i] += yi[k] * cc[i] * w;
		    }
	    }

	/*--------- Solve Ax=b ----------*/
	{
	    if (!isfinite(A(0,0)))
		{
		    //            DebugBreak();
		    return false;
		}
	    else
		{
		    gtb::AVec<double> x(NCF);
	           
		    SVDSolve(A, b, x);

		    /*
		     * Copy result to the coeeficients of the polynom
		     */
		    coefs[0] = x[2];
		    coefs[1] = x[1];
		    coefs[2] = x[0];
		}
	    return true;
	}
    }



    real_type coefs[3];
};


real_type edge_curvature(const vector<Point3> &points, const vector<real_type> &weights, int si) {

    vector<real_type> ss;
    vector<real_type> xs, ys, zs;
    for (unsigned i=0; i<points.size(); i++) {

	if (i==0)
	    ss.push_back(0);
	else
	    ss.push_back(ss.back() + Point3::distance(points[i-1], points[i]));

	xs.push_back(points[i][0]);
	ys.push_back(points[i][1]);
	zs.push_back(points[i][2]);
    }

    for (unsigned i=0; i<ss.size(); i++) {
	ss[i] /= ss.back();
    }


    Quadratic pcurves[3];

    if (!pcurves[0].fit(ss, xs, weights)) return 0;
    if (!pcurves[1].fit(ss, ys, weights)) return 0;
    if (!pcurves[2].fit(ss, zs, weights)) return 0;


    real_type s = ss[si];

    Vector3 qd(pcurves[0].eval1(s), pcurves[1].eval1(s), pcurves[2].eval1(s));
    Vector3 qdd(pcurves[0].eval2(s), pcurves[1].eval2(s), pcurves[2].eval2(s));

    real_type k = qd.cross(qdd).length() / gtb::ipow(qd.length(), 3);
    return k;
}



///////////////////////////////////////////////////////////////////////////////
// MeshProjector
MeshProjector::kdtree_type *MeshProjector::GetKdTree()
{
    return kdtree;
}

MeshProjector::MeshProjector(const TriangleMesh &m, const vector<int> *pointsides) : mesh(m) {

    vector<int> fps;
    for (unsigned f=0; f<mesh.faces.size(); f++) {
	// for csg operations, we can skip a bunch of the faces that never see the light of day
	if (!pointsides ||
	    (*pointsides)[m.faces[f].verts[0]]>0 || 
	    (*pointsides)[m.faces[f].verts[1]]>0 || 
	    (*pointsides)[m.faces[f].verts[2]]>0)
	    fps.push_back(f);
    }

    kdtree = new kdtree_type(fps, *this);
}

MeshProjector::~MeshProjector() {
    delete kdtree;
}



int MeshProjector::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const {
	
    kdtree_type::OrderedTraverse ot(*kdtree, fp, *this);
    int fi; ot.next(fi);
    const TriangleMeshFace &f = mesh.faces[fi];
    Triangle3 tri(mesh.verts[f.verts[0]].point, mesh.verts[f.verts[1]].point, mesh.verts[f.verts[2]].point);
    tp = tri.closest_point(fp);

    real_type b[3];
    tri.get_barycentric_coordinates(tp, b[0], b[1], b[2]);

    tn = Vector3(0,0,0);
    tn += b[0] * mesh.verts[f.verts[0]].normal;
    tn += b[1] * mesh.verts[f.verts[1]].normal;
    tn += b[2] * mesh.verts[f.verts[2]].normal;
    tn.normalize();

    return PROJECT_SUCCESS;
}



int MeshProjector::ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const {

#ifdef CLOSEST_POINT_PROJECTION
    return ProjectPointClosest(fp, tp, tn);
#else


    Vector3 axis = base2.position - base1.position;
    axis.normalize();

    Point3 center = base1.position + (axis.dot(fp - base1.position)) * axis;
    real_type rad = Point3::distance(center, fp);


    Vector3 extremes;
    for (int i=0; i<3; i++) {
		extremes[i] = fabs(rad * cos(M_PI/2 - acos(axis[i]))) + 1e-5;
    }

    Box3 ringbox(center-extremes, center+extremes);


    static_vector(int,possible_intersects);
    kdtree->GetIntersectedBoxes(*this, ringbox, possible_intersects);


    real_type closest_intersect = 1e34;

    for (unsigned i=0; i<possible_intersects.size(); i++) {

		const TriangleMeshFace &f = mesh.faces[possible_intersects[i]];
		Triangle3 tri(mesh.verts[f.verts[0]].point, mesh.verts[f.verts[1]].point, mesh.verts[f.verts[2]].point);

		if (tri.area() == 0) continue;

		Vector3 fnorm = tri.normal();

		real_type d = fnorm.dot(center) - fnorm.dot(tri.A());

		if (fabs(d) > rad) continue;

		Vector3 v = fnorm.cross(axis);
		v.normalize();

		real_type x = sqrt(rad*rad - d*d);

		for (int j=0; j<2; j++) {

			Point3 p = center + -d*fnorm + ((j==0)?-x:x) * v;

			// see if it's actually inside the triangle
			bool inside = true;
			for (int k=0; k<3; k++) {
				if (fnorm.dot((tri[(k+1)%3]-tri[k]).cross(p-tri[k])) < 0) {
					inside = false; break;
				}
			}
			if (!inside) continue;

			if (axis.dot(fnorm.cross(p-center)) > 0)
				continue;

			real_type tdist = Point3::distance(p, fp);
			if (tdist < closest_intersect) {
				closest_intersect = tdist;
				tp = p;

				real_type b[3];
				tri.get_barycentric_coordinates(tp, b[0], b[1], b[2]);

				tn = Vector3(0,0,0);
				tn += b[0] * mesh.verts[f.verts[0]].normal;
				tn += b[1] * mesh.verts[f.verts[1]].normal;
				tn += b[2] * mesh.verts[f.verts[2]].normal;
				tn.normalize();
			}

		}

    }

	if (closest_intersect < 1e30)
		return PROJECT_SUCCESS;
	return PROJECT_FAILURE;

#endif

}


Box3 MeshProjector::bounding_box(int face) const {
    const TriangleMeshFace &f = mesh.faces[face];
    return Box3::bounding_box(mesh.verts[f.verts[0]].point, mesh.verts[f.verts[1]].point, mesh.verts[f.verts[2]].point);
}

real_type MeshProjector::distance(int face, const Point3 &from) const {
    const TriangleMeshFace &f = mesh.faces[face];
    Triangle3 tri(mesh.verts[f.verts[0]].point, mesh.verts[f.verts[1]].point, mesh.verts[f.verts[2]].point);
    return Point3::distance(from, tri.closest_point(from));
}



// try to find the boundaries in a triangle soup
void MeshProjector::FindSoupBoundaries(vector< vector<int> > &boundaries, real_type d) const {

    // get candidate edges:
    // edges that when a point is projected from off to it's side
    // get projected back onto the edge
    std::set< std::pair<int,int> > bedges;
    for (unsigned t=0; t<mesh.faces.size(); t++) {
	for (int i=0; i<3; i++) {
	    if (mesh.faces[t].nbrs[i] < 0) {
		// try projecting a point out to the side onto the mesh
		Point3 v1 = mesh.verts[mesh.faces[t].verts[i]].point;
		Point3 v2 = mesh.verts[mesh.faces[t].verts[(i+1)%3]].point;
		Vector3 n1 = mesh.verts[mesh.faces[t].verts[i]].normal;
		Vector3 n2 = mesh.verts[mesh.faces[t].verts[(i+1)%3]].normal;
		Vector3 norm = (n1+n2).normalized();
		Vector3 perp = (v2-v1).cross(norm).normalized();
		Point3 mid = v1 + (real_type)0.5 * (v2-v1);

		real_type td = d;
		if (td < 0) {
		    td = -td*Point3::distance(v1,v2);
		}

		Point3 fp = mid + td*perp;
		Point3 tp;
		Vector3 tn;

		if (ProjectPoint(fp, tp, tn) == PROJECT_SUCCESS) {
		    if (Point3::distance(mid, tp) < td*0.1) {
			bedges.insert(std::pair<int,int>(mesh.faces[t].verts[i], mesh.faces[t].verts[(i+1)%3]));
		    }
		}
	    }
	}
    }
#if 0
    dbgClear();
    for (std::set<std::pair<int,int> >::iterator i=bedges.begin();
	 i != bedges.end(); ++i) {

	vector<Point3> line(2);
	line[0] = mesh.verts[i->first].point;
	line[1] = mesh.verts[i->second].point;
	DbgPLines::add(line, 0);
    }
    redrawAndWait(' ');
#endif


    // try stitching them together into loops
    while (1) {

	if (bedges.empty())
	    break;


	vector<int> loop;
	loop.push_back(bedges.begin()->first);
	loop.push_back(bedges.begin()->second);
	bedges.erase(bedges.begin());

	bool made_loop = false;
	while (1) {

	    if (bedges.empty())
		break;

	    // look for the edge at the end of the loop
	    real_type dist = 1e34;
	    std::set<std::pair<int,int> >::iterator best = bedges.end();

	    for (std::set<std::pair<int,int> >::iterator i=bedges.begin();
		 i != bedges.end(); ++i) {

		real_type tdist = Point3::distance(mesh.verts[loop.back()].point, mesh.verts[i->first].point);
		if (tdist < dist) {
		    dist = tdist;
		    best = i;
		}
	    }

	    if (best == bedges.end()) {
		cerr<<"no best edge??"<<endl;
		break;
	    }

	    real_type elen = Point3::distance(mesh.verts[best->first].point, mesh.verts[best->second].point);

	    // skipped to an edge too far away
	    if (dist > elen*0.1) {
		cerr<<"dist too far from last point in loop"<<endl;
		break;
	    }


	    // see if we closed the loop
	    real_type clen = Point3::distance(mesh.verts[loop.front()].point, mesh.verts[best->second].point);
	    if (clen < elen*0.1) {
		cerr<<"loop done"<<endl;
		made_loop = true;
		break;
	    }

#if 0
	    vector<Point3> line(2);
	    line[0] = mesh.verts[best->first].point;
	    line[1] = mesh.verts[best->second].point;
	    DbgPLines::add(line, 0);
	    redrawAndWait(' ');
#endif

	    // add to the loop
	    loop.push_back(best->second);
	    bedges.erase(best);
	}


	if (made_loop) {
	    boundaries.push_back(loop);
#if 0
	    dbgClear();
	    vector<Point3> line;
	    for (unsigned i=0; i<loop.size(); i++) {
		line.push_back(mesh.verts[loop[i]].point);
	    }
	    line.push_back(mesh.verts[loop[0]].point);
	    DbgPLines::add(line, 1);
	    redrawAndWait(' ');
#endif
	}
    }
}




///////////////////////////////////////////////////////////////////////////////
// MeshGuidanceField

MeshGuidanceField::kdtree_type *MeshGuidanceField::GetKdTree()
{
    return kdtree;
}

gtb::tPoint3<double> toDouble(const Point3 &p) {
    return gtb::tPoint3<double>(p[0], p[1], p[2]);
}
gtb::tVector3<double> toDouble(const Vector3 &v) {
    return gtb::tVector3<double>(v[0], v[1], v[2]);
}


real_type SubdivCurvature(const Point3 &_center, const vector<Point3> &_ring, bool edge, int niter) {

    assert(niter>0);

    const unsigned N = _ring.size();
    gtb::tPoint3<double> center[2];
    vector< gtb::tPoint3<double> > rings[2];

    rings[0].resize(N);
    rings[1].resize(N);

    int s=0;
    center[0] = toDouble(_center);
    for (unsigned i=0; i<N; i++) {
	rings[0][i] = toDouble(_ring[i]) - center[0];
    }
    center[0] = gtb::tPoint3<double>(0,0,0);


    double beta = (N>3) ? (3.0/(8.0*N)) : (3.0/16.0);


    for (int iter=0; iter<niter; iter++) {

	int os=1-s;

	center[os] = gtb::tPoint3<double>(0,0,0);
	if (edge) {
	    center[os].add_scaled(center[s], 6.0/8.0);
	    center[os].add_scaled(rings[s].front(), 1.0/8.0);
	    center[os].add_scaled(rings[s].back(), 1.0/8.0);
	} else {
	    center[os].add_scaled(center[s], 1.0-N*beta);
	    for (unsigned i=0; i<N; i++) {
		center[os].add_scaled(rings[s][i], beta);
	    }
	}

	for (unsigned i=0; i<N; i++) {
	    rings[os][i] = gtb::tPoint3<double>(0,0,0);
	    if (edge && (i==0 || i==N-1)) {
		rings[os][i].add_scaled(center[s], 4.0/8.0);
		rings[os][i].add_scaled(rings[s][i], 4.0/8.0);
	    } else {
		rings[os][i].add_scaled(center[s], 3.0/8.0);
		rings[os][i].add_scaled(rings[s][i], 3.0/8.0);
		rings[os][i].add_scaled(rings[s][(i+1)%N], 1.0/8.0);
		rings[os][i].add_scaled(rings[s][(i+N-1)%N], 1.0/8.0);
	    }
	}

	s = os;
    }


    // compute the normal
    gtb::tVector3<double> normal(0,0,0);
    for (unsigned i=0; i<((edge)?N-1:N); i++) {
	normal += (rings[s][i] - center[s]).cross(rings[s][(i+1)%N] - center[s]).normalized();
    }
    normal.normalize();


    // stuff it into this surfelset crap so we can fit a polynomial to it
    gtb::tsurfel_set<double> nbhdset;
    gtb::tsurfelset_view<double> nbhdview(nbhdset);
    gtb::tPlane<double> plane(normal, toDouble(Point3(0,0,0)));
    gtb::plane_transformation<double> T(plane, center[s]);
    for (unsigned i=0; i<N; i++) {
	nbhdview.insert(nbhdset.size());
	nbhdset.insert_vertex(rings[0][i]);
	nbhdview.insert(nbhdset.size());
	nbhdset.insert_vertex(rings[1][i]);
    }

    mls::Poly2<double> poly;
    gtb::tsurfel_set<double> stdpts;
    mls::GaussianWeightFunction<double> wf;	wf.set_radius(1);

    mls::CProjection<double>::WeightedPolyFit(nbhdview, T, center[s], &wf, &poly, stdpts);

    double k1, k2;
    poly.curvatures(0,0,k1,k2);

	return k2;
}



void SolveNormalEquations(const AMat &nA, const AVec &nB, gtb::AVec<double> &x) {

    gtb::AMat<double> A(nA.columns(), nA.columns());

    for (int m=0; m<A.rows(); m++) {
	for (int n=0; n<A.columns(); n++) {
	    A(m,n) = 0;
	    for (int i=0; i<nA.rows(); i++) {
		A(m,n) += (double)nA(i,m)*(double)nA(i,n);
	    }
	}
    }

    gtb::AVec<double> b(nA.columns());
    for (int m=0; m<b.size(); m++) {
	b[m] = 0;
	for (int i=0; i<nA.rows(); i++) {
	    b[m] += (double)nB[i] * (double)nA(i,m);
	}
    }

    x.resize(b.size());
    SVDSolve(A, b, x);
}


real_type ExtendedQuadricCurvature(const Point3 &center, const vector<Point3> &points, const Vector3 &first_normal) {

    Vector3 normal, lnormal;
    normal = first_normal;

    AMat A(points.size(), 5);
    AVec b(points.size());
    gtb::AVec<double> x(5);

    int iter=0;
    do {
	iter++;

	Plane plane(normal, Point3(0,0,0));
	plane_transformation T(plane, center);

	Point3 _center = T.ToPlane(center);

	for (unsigned p=0; p<points.size(); p++) {
	    Point3 local = T.ToPlane(points[p]);
	    A(p,0) = local[0]*local[0];
	    A(p,1) = local[0]*local[1];
	    A(p,2) = local[1]*local[1];
	    A(p,3) = local[0];
	    A(p,4) = local[1];
	    b[(int)p] = local[2];
	}

	SolveNormalEquations(A, b, x);

	lnormal = normal;
	normal = T.PlaneNormal2World(Vector3(-x[3],-x[4],1).normalized());

    } while (iter<10 && normal.dot(lnormal) < 0.99);



    double descrim = std::max(0.0, (x[0]-x[2])*(x[0]-x[2]) + x[1]*x[1]);

	if (x[0]+x[2] > 0)
	    return (real_type)(x[0] + x[2] + sqrt(descrim));
	else
	    return (real_type)(x[0] + x[2] - sqrt(descrim));
}







MeshGuidanceField::~MeshGuidanceField() {
    delete kdtree;
}

MeshGuidanceField::MeshGuidanceField(int curv_sub, const TriangleMesh &mesh, real_type rho, real_type min_step, real_type max_step, real_type reduction)
    : GuidanceField(rho, min_step, max_step, reduction), kdGetPoint(mesh) {


    // setup the kdtree
    kdtree = new kdtree_type(10, mesh.bounding_box(), kdGetPoint);
    for (unsigned i=0; i<mesh.verts.size(); i++)
        kdtree->Insert(i, false);
    kdtree->MakeTree();


    mls::GaussianWeightFunction<double> wf;
    wf.set_radius(1);


    TriangleMesh cmesh = mesh;
    for (int i=0; i<curv_sub; i++) {
	cmesh.LoopSubdivide();
    }


    for (unsigned v=0; v<mesh.verts.size(); v++) {


	// fit a polynomial
	const int numrings = 3;
	vector< std::set<int> > rings(numrings+1);
	rings[0].insert(v);

	for (int r=0; r<numrings; r++) {
	    rings[r+1] = rings[r];

	    for (std::set<int>::iterator i=rings[r].begin(); i!=rings[r].end(); ++i) {

		for (TriangleMesh::VertexVertexIteratorI vi(cmesh, *i); !vi.done(); ++vi) {
		    rings[r+1].insert(*vi);
		}
	    }
	}

	vector<Point3> nbrhood;
	Vector3 normest = cmesh.verts[v].normal;
	for (unsigned r=0; r<rings.size(); r++) {
	    for (std::set<int>::iterator i=rings[r].begin(); i!=rings[r].end(); ++i) {
		nbrhood.push_back(cmesh.verts[*i].point);

		if (cmesh.verts[*i].normal.length() > 0.5)
			normest += cmesh.verts[*i].normal;
	    }
	}
	normest.normalize();


	real_type kmax;
	if (normest.length() > 0.5)
	    kmax = ExtendedQuadricCurvature(cmesh.verts[v].point, nbrhood, normest);
	else 
		kmax = 0;


	if (cmesh.VertOnBoundary(v)) {

	    const int num =2;

	    vector<int> eptsi;
	    eptsi.push_back(v);
            
	    for (int i=0; i<num; i++) {
		TriangleMesh::VertexVertexIteratorI vi(cmesh, eptsi.back());
		eptsi.push_back(*vi);
	    }

	    reverse(eptsi);


	    for (int i=0; i<num; i++) {
		TriangleMesh::VertexVertexIteratorI vi(cmesh, eptsi.back());
		eptsi.push_back(0);
		while (!vi.done()) {
		    eptsi.back() = *vi;
		    ++vi;
		}
	    }


	    vector<Point3> epts;
	    vector<real_type> weights;
	    for (unsigned i=0; i<eptsi.size(); i++) {
		epts.push_back(cmesh.verts[eptsi[i]].point);
		weights.push_back(1);
	    }


	    double k = edge_curvature(epts, weights, epts.size()/2);

	    if (fabs(k) > fabs(kmax))
		kmax = k;
	}


	ideal_length.push_back(MaxCurvatureToIdeal(kmax));



	if ((int)(v * 100.0 / mesh.verts.size()) != (int)((v+1) * 100.0 / mesh.verts.size())) {
	    cerr<<"                  \r"<<((int)((v+1) * 100.0 / mesh.verts.size()))<<"%"<<std::flush;
	}

    }
    cerr << std::endl;


	int num_orig = ideal_length.size();

    vector<int> marked(ideal_length.size()); // for removal
    Trimmer<GetPointMesh> trimmer(*this, kdGetPoint, ideal_length);
    trimmer.Trim(marked);


    // rebuild the kdtree
    delete kdtree;

	int num_trimmed = 0;
    kdtree = new kdtree_type(10, mesh.bounding_box(), kdGetPoint);
    for (unsigned i=0; i<mesh.verts.size(); i++) {
		if (!marked[i]) {
			kdtree->Insert(i, false);
			num_trimmed++;
		}
    }
    kdtree->MakeTree();

	
	cerr<<num_orig<<" original guidance field points"<<endl;
	cerr<<num_trimmed<<" guidance field points after trimming"<<endl;
}



void* MeshGuidanceField::OrderedPointTraverseStart(const Point3 &p) {

	//    kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = 
	return new kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
	
}

bool MeshGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {

	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = kdGetPoint(n);
    ideal = ideal_length[n];
    return true;
}

void MeshGuidanceField::OrderedPointTraverseEnd(void *ctx) {
	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
    delete kdOrderedTraverse;
}


void MeshGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad) {
    pts.resize(kdGetPoint.mesh.verts.size());
    norms.resize(kdGetPoint.mesh.verts.size());
    rad.resize(kdGetPoint.mesh.verts.size());
    for (unsigned i=0; i<kdGetPoint.mesh.verts.size(); i++) {
	pts[i]   = kdGetPoint.mesh.verts[i].point;
	norms[i] = kdGetPoint.mesh.verts[i].normal;
	rad[i]   = 1 / ideal_length[i];
    }
}


