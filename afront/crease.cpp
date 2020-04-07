
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
#include "crease.h"

using namespace std;

//#define DEBUG_CREASES

void CreaseExtractor::MakeGraph()
{
    total_half_edges = 0;

    // First, threshold all creases
    for (unsigned fi=0; fi<m.faces.size(); ++fi) {
	const gtb::TriangleMeshFace &f = m.faces[fi];
	Vector3 fn = m.FaceNormal(f);
	for (int ei=0; ei<3; ++ei) {
	    if (f.nbrs[ei] == -1)
		continue;
	    int fi2 = f.nbrs[ei];
	    const gtb::TriangleMeshFace &nbr = m.faces[fi2];
	    Vector3 nbrn = m.FaceNormal(nbr);
	    if (fn.dot(nbrn) < sharp) {
		//		cerr << f.verts[ei] << " -> " << f.verts[(ei+1)%3] << endl;
		creases[f.verts[ei]].push_back(f.verts[(ei+1)%3]);
		total_half_edges++;
	    }
	}
    }

    // Now we find the kink vertices
    map<int, list<int> >::iterator 
	b = creases.begin(),
	e = creases.end();

    for (; b!=e; ++b)
	if (b->second.size() > 2)
	    kink_vertices.push_back(b->first);
}

//! \brief Returns a vertex from the graph that's also a kink vertex
int CreaseExtractor::FindKinkVertex(const Graph &c)
{
    list<int>::iterator 
	b = kink_vertices.begin(),
	e = kink_vertices.end();
    while (b != e) {
	if (c.find(*b) != c.end())
	    return *b;
	++b;
    }
    return -1;
}

bool CreaseExtractor::IsKinkVertex(int v)
{
    return find(kink_vertices.begin(),
		kink_vertices.end(),
		v) != kink_vertices.end();
}

#define CIRCULAR_NEXT(itor, mesh, vertex)				\
{									\
    ++itor;								\
    if (itor.done())							\
	itor = TriangleMesh::VertexVertexIteratorI(mesh,vertex);	\
}


void CreaseExtractor::ComputeKinks()
{
    Graph creases_copy = creases;
    int processed_half_edges = 0;
    while (processed_half_edges < total_half_edges) {
	// try to start at a kink vertex
	int start = FindKinkVertex(creases_copy);
	if (start == -1) {
	    // There was no kink vertex, choose an arbitrary one 
	    // (it will necessarily be in the middle of a kink loop)
	    start = creases_copy.begin()->first;
	    // we arbitrarily name this vertex a kink vertex
	    kink_vertices.push_back(start);
	}

	int current = start;
	vector<int> kink;
	kink.push_back(start);
	list<int>::iterator itor = creases_copy[current].begin();
	int next = *itor;
	// Traverse the kink all the way to a vertex
	while (!IsKinkVertex(next)) {
	    kink.push_back(next);
	    TriangleMesh::VertexVertexIteratorI vvi(m, next);
	    while (*vvi != current) ++vvi;
	    list<int>::iterator f;
	    // finds the next crease vertex on the star of the current vertex
	    do {
		CIRCULAR_NEXT(vvi, m, next);
		list<int>::iterator
		    e = creases_copy[next].end();
		f = find(creases_copy[next].begin(), e, *vvi);
		if (f != e)
		    break;
	    } while (1);
	    current = next;
	    next = *vvi;
	}
	kink.push_back(next);
#ifdef DEBUG_CREASES
	{
	    vector<Point3> tv(kink.size());
	    for (unsigned i=0; i<kink.size(); ++i)
		tv[i] = m.verts[kink[i]].point;
 	    DbgPLines::add(tv);
 	    redrawAndWait(' ');
	}
#endif
	// Erase crease half-edges that were put in the kink
	for (unsigned i=0; i<kink.size()-1; ++i) {
	    creases_copy[kink[i]].erase(find(creases_copy[kink[i]].begin(),
					     creases_copy[kink[i]].end(),
					     kink[i+1]));
	    processed_half_edges++;
 	    if (!creases_copy[kink[i]].size())
		creases_copy.erase(kink[i]);
	}
	// If kink is a twig, artificially split it at the end of the twig
	set<int> vs;
	bool is_twig = false;
	int twig_vertex = -1;
	for (unsigned i=0; i<kink.size(); ++i) {
	    if ((vs.find(kink[i]) != vs.end()) && // kink meets itself
		(i < kink.size()-1)) {            // but it's not a loop kink
		is_twig = true;
		twig_vertex = i-1;
		break;
	    }
	    vs.insert(kink[i]);
	}
	if (is_twig) {
	    vector<int> twig_to_end, twig_from_end;
	    for (int i=0; i<=twig_vertex; ++i)
		twig_to_end.push_back(kink[i]);
	    for (unsigned i=twig_vertex; i<kink.size(); ++i)
		twig_from_end.push_back(kink[i]);
	    kinks.push_back(twig_to_end);
	    kinks.push_back(twig_from_end);
	} else
	    kinks.push_back(kink);
    }
}

Vector3 CreaseExtractor::ComputeEdgeNormal(vector<int> &kink, int index)
{
    const Point3 &indexPoint = m.verts[kink[index]].point;
    Point3 prevPoint = m.verts[kink[index-1]].point;

    TriangleMesh::VertexVertexIteratorI vvi(m, kink[index]);
    while (*vvi != kink[index-1])
	CIRCULAR_NEXT(vvi, m, kink[index]);
    CIRCULAR_NEXT(vvi, m, kink[index]);
    Point3 nextPoint = m.verts[*vvi].point;

    Vector3 accNormal(0,0,0);
    Vector3 da = prevPoint - indexPoint,
	db = nextPoint - indexPoint,
	normal = da.cross(db);
    normal.normalize();
    accNormal += normal;

    do {
	prevPoint = nextPoint;
	CIRCULAR_NEXT(vvi, m, kink[index]);
	nextPoint = m.verts[*vvi].point;
	da = prevPoint - indexPoint;
	db = nextPoint - indexPoint;
	normal = da.cross(db);
	normal.normalize();
	accNormal += normal;
    } while (*vvi != kink[index+1]);
    accNormal.normalize();
#ifdef DEBUG_CREASES
    {
	vector<Point3> p;
	p.push_back(m.verts[kink[index]].point);
	p.push_back(m.verts[kink[index]].point +
		    (accNormal * (m.verts[kink[index]].point - 
				  m.verts[kink[index-1]].point).length()));
	DbgPLines::add(p);
    }
#endif
    return accNormal;
}

void CreaseExtractor::ComputeResampledKinks()
{
    for (unsigned i=0; i<kinks.size(); ++i) {
	pair<int,int> key(min(kinks[i].front(), kinks[i].back()),
			  max(kinks[i].front(), kinks[i].back()));
	ResampledCurve value;
	map<pair<int, int>, ResampledCurve>::iterator 
	    f = resampledKinks.find(key),
	    e = resampledKinks.end();
	if (f == e) {
	    set<int> points;
	    vector<Point3> ploop;
	    vector<vector<Vector3> > nloop(1);
	    if (key.first == kinks[i].front()) {
		// Traverse in regular order
		for (unsigned j=0; j<kinks[i].size(); ++j) {
		    ploop.push_back(m.verts[kinks[i][j]].point);
		    nloop[0].push_back(m.verts[kinks[i][j]].normal);

//  		    if (kinks[i].size() == 2)
// 			nloop[0].push_back(m.verts[kinks[i][j]].normal);
// 		    else {
// 			unsigned v = j;
// 			v = max((unsigned)1, v);
// 			v = min(kinks[i].size()-2, v);
// 			nloop[0].push_back(ComputeEdgeNormal(kinks[i], v));
// 		    }
		}
	    } else {
		// Traverse in reverse order
		for (int j=kinks[i].size()-1; j>=0; --j) {
		    ploop.push_back(m.verts[kinks[i][j]].point);
		    nloop[0].push_back(m.verts[kinks[i][j]].normal);
		    
// 		    if (kinks[i].size() == 2)
// 			nloop[0].push_back(m.verts[kinks[i][j]].normal);
// 		    else {
// 			unsigned v = j;
// 			v = max((unsigned)1, v);
// 			v = min(kinks[i].size()-2, v);
// 			nloop[0].push_back(ComputeEdgeNormal(kinks[i], v));
// 		    }
		}
	    }
	    vector<vector<Vector3> > nlooprs(1);
	    guidance->ResampleCurve(ploop, nloop, value.points, nlooprs, true);
	    value.normals = nlooprs[0];
	    resampledKinks[key] = value;
	}
    }

    // Go over resampled kinks and compute the vertex indices 
    // at the same time, write on crease_points and crease_onormals

    map<pair<int,int>, ResampledCurve>::iterator
	b = resampledKinks.begin(),
	e = resampledKinks.end();
    for (; b!=e; ++b) {
	for (unsigned i=0; i<b->second.points.size(); ++i) {
	    b->second.indices.push_back(crease_points.size());
	    crease_points.push_back(b->second.points[i]);
	    crease_onormals.push_back(b->second.normals[i]);
	}
    }
}

int CreaseExtractor::FindKinkWithEndPoint(int v, int going_through)
{
    for (unsigned i=0; i<kinks.size(); ++i) {
	if (v == kinks[i].front() &&
	    (going_through == -1 ||
	     going_through == kinks[i][1]))
	    return i;
	if (v == kinks[i].back() &&
	    (going_through == -1 ||
	     going_through == kinks[i][kinks[i].size()-2]))
	    return i;
    }
    return -1;    
}

int CreaseExtractor::OtherKinkEndPoint(int kink_index, int endpoint)
{
    if (kinks[kink_index].front() == endpoint)
	return kinks[kink_index].back();
    else
	return kinks[kink_index].front();
}

int CreaseExtractor::KinkVertexNeighbor(int kink_index, int endpoint)
{
    if (kinks[kink_index].front() == endpoint)
	return kinks[kink_index][1];
    else
	return kinks[kink_index][kinks[kink_index].size()-2];
}

void CreaseExtractor::ReverseFrontIfNecessary
(vector<int> &front,
 vector<Vector3> &normals)
{
    vector<int> result_front;
    vector<Vector3> result_normals;

    copy(front.rbegin(), front.rend(), back_inserter(result_front));
    copy(normals.rbegin(), normals.rend(), back_inserter(result_normals));

    normals = result_normals;
    front = result_front;
}

void CreaseExtractor::AddKinkToFront
(vector<int> &front, vector<Vector3> &normals,
 int kink_index, int endpoint)
{
    int from = endpoint,
	to = OtherKinkEndPoint(kink_index, endpoint);
    pair<int,int> resampled_key(min(from, to),
				max(from, to));
    ResampledCurve &c = resampledKinks[resampled_key];
    if (from == min(from, to)) {
#ifdef DEBUG_CREASES
 	for (unsigned i=0; i<c.points.size()-1; ++i)
 	    DbgPoints::add(c.points[i], 1,0,0);
 	redrawAndWait(' ');
#endif
	copy(c.indices.begin(), c.indices.end()-1, back_inserter(front));
	copy(c.normals.begin(), c.normals.end()-1, back_inserter(normals));
    } else {
#ifdef DEBUG_CREASES
 	for (unsigned i=1; i<c.points.size(); ++i)
 	    DbgPoints::add(c.points[i], 1,0,0);
 	redrawAndWait(' ');
#endif
	copy(c.indices.rbegin(), c.indices.rend()-1, back_inserter(front));
	copy(c.normals.rbegin(), c.normals.rend()-1, back_inserter(normals));
    }
}

void CreaseExtractor::ComputeFronts()
{
    Graph creases_copy = creases;
    int processed_half_edges = 0;
    map<int, bool> loop_kink_reversal;
    while (processed_half_edges < total_half_edges) {
	int start = FindKinkVertex(creases_copy);
	vector<int> front;
	vector<Vector3> normals;
	vector<int> kinks_used;
	vector<bool> kinks_reversed;
	int current = start;
	int going_through = *creases_copy[current].begin();
	int kink_index = FindKinkWithEndPoint(current, going_through);
	do {
	    AddKinkToFront(front, normals, kink_index, current);
	    kinks_used.push_back(kink_index);
	    kinks_reversed.push_back(current == kinks[kink_index].back());
	    int next = OtherKinkEndPoint(kink_index, current);
	    if (next == start)
		break;
	    TriangleMesh::VertexVertexIteratorI vvi(m, next);
	    int pred = KinkVertexNeighbor(kink_index, next);
	    while (*vvi != pred) ++vvi;
	    list<int>::iterator f;
	    // finds the next crease vertex on the star of the current vertex
	    do {
		CIRCULAR_NEXT(vvi, m, next);
		list<int>::iterator
		    e = creases_copy[next].end();
		f = find(creases_copy[next].begin(), e, *vvi);
		if (f != e)
		    break;
	    } while (1);
	    current = next;
	    kink_index = FindKinkWithEndPoint(current, *vvi);
	} while (1);
	if (kinks_used.size() == 1) {
	    // Special ugly case for loop kink
	    if (loop_kink_reversal[kinks_used[0]] == false) {
		loop_kink_reversal[kinks_used[0]] = true;
		kinks_reversed[0] = false;
		ReverseFrontIfNecessary(front, normals);
	    } else
		kinks_reversed[0] = true;
	} else
	    ReverseFrontIfNecessary(front, normals);

	for (unsigned i=0; i<kinks_used.size(); ++i) {
	    vector<int> &kink = kinks[kinks_used[i]];
	    if (!kinks_reversed[i]) {
		for (unsigned j=0; j<kink.size()-1; ++j) {
		    creases_copy[kink[j]].
			erase(find(creases_copy[kink[j]].begin(),
				   creases_copy[kink[j]].end(),
				   kink[j+1]));
		    processed_half_edges++;
		    if (!creases_copy[kink[j]].size())
			creases_copy.erase(kink[j]);
		}
	    } else {
		for (unsigned j=kink.size()-1; j>=1; --j) {
		    creases_copy[kink[j]].
			erase(find(creases_copy[kink[j]].begin(),
				   creases_copy[kink[j]].end(),
				   kink[j-1]));
		    processed_half_edges++;
		    if (!creases_copy[kink[j]].size())
			creases_copy.erase(kink[j]);
		}
	    }
	}
	kinks_used.clear();
	kinks_reversed.clear();

	crease_indices.push_back(front);
	crease_normals.push_back(normals);
    }
}


void CreaseExtractor::AddCornerTriangles()
{
    MeshProjector::kdtree_type *kd = projector->GetKdTree();
    for (unsigned i=0; i<crease_indices.size(); ++i) {
	std::vector<int> &front = crease_indices[i];
// 	bool cut_vertices;
// 	do {
// 	    cut_vertices = false;
	    for (unsigned j=0; j<front.size(); ++j) {
		int v1 = j,
		    v2 = (j+1) % front.size(),
		    v3 = (j+2) % front.size();
		Vector3 
		    e1 = crease_points[front[v2]] - crease_points[front[v1]],
		    e2 = crease_points[front[v3]] - crease_points[front[v2]];

		float l = e1.length();
		e1.normalize();
		e2.normalize();



		if (e1.dot(e2) >= 0.2)
		    continue;
		{
		    vector<Point3> v;
		    v.push_back(crease_points[front[v1]]);
		    v.push_back(crease_points[front[v2]]);
		    v.push_back(crease_points[front[v3]]);
		    DbgPLines::add(v);
		    DbgPoints::add(crease_points[front[v1]], 1, 0, 0);
		    DbgPoints::add(crease_points[front[v2]], 1, 0, 0);
		    DbgPoints::add(crease_points[front[v3]], 1, 0, 0);
		    //		    cerr << e1.dot(e2) << endl;
		    //		    redrawAndWait(' ');
		}
		//		cerr << "Found candidate" << endl;
		Vector3 lv(l,l,l);
		Vector3 cp = e1.cross(e2);
		
	    
		Box3 box(crease_points[front[v2]]-lv,
			 crease_points[front[v2]]+lv);
		static_vector(int,possible_intersects);
		kd->GetIntersectedBoxes(*projector, box, possible_intersects);
		bool right_side = false;
		//		cerr << possible_intersects.size() << endl;
		for (unsigned fi=0; fi<possible_intersects.size(); ++fi) {
		    const TriangleMeshFace &f = m.faces[possible_intersects[fi]];
		    Triangle3 tri(m.verts[f.verts[0]].point,
				  m.verts[f.verts[1]].point,
				  m.verts[f.verts[2]].point);
		    Vector3 fnorm = tri.normal();
		    real_type d = fnorm.dot(cp);
		    //		    cerr << cp << " - " << fnorm << " - " << d << endl;
		    if (d > 0.5) {
			right_side = true;
			break;
		    }
		}

		if (right_side) {
		    //		    cerr << "will cut ear" << endl;
		    // Cut ear
		    corner_triangles.push_back(front[v1]);
		    corner_triangles.push_back(front[v2]);
		    corner_triangles.push_back(front[v3]);
		    //		    cut_vertices = true;
		    front.erase(front.begin() + v2);
		    j--;
		}
	    }
// 	} while (cut_vertices);
    }
}

void CreaseExtractor::DeleteSmallCreases()
{
    for (unsigned i=0; i<crease_indices.size(); ++i) {
	cerr << crease_indices[i].size() << endl;
	if (crease_indices[i].size() < small_crease_threshold) {
	    crease_indices.erase(crease_indices.begin() + i);
	    crease_normals.erase(crease_normals.begin() + i);
	    --i;
	}
    }
}

void CreaseExtractor::Extract()
{
    cerr << "Computing creases.." << endl;
    MakeGraph();
    if (!creases.size())
	return;
    cerr << "Computing kinks.." << endl;
    ComputeKinks();
    cerr << "Resampling kinks.." << endl;
    ComputeResampledKinks();
    cerr << "Computing front.." << endl;
    ComputeFronts();
    cerr << "Deleting small creases..." << endl;
    DeleteSmallCreases();
    cerr << "Adding corner triangles..." << endl;
    AddCornerTriangles();
   
    cerr << "Ok." << endl;
}
