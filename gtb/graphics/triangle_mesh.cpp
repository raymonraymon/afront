
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


#include <gtb/gtb.hpp>
#include <gtb/graphics/triangle_mesh.h>

using std::vector;
using std::cerr;
using std::endl;

GTB_BEGIN_NAMESPACE

template <class T>
tTriangleMesh<T>::tTriangleMesh() 
{
	tModel<T>::invalidate_all();
}



template <class T>
tTriangleMesh<T>::~tTriangleMesh() 
{
}


// read a Hugues Hoppe .m format mesh file
template <class T>
bool tTriangleMesh<T>::ReadHHM(const char *fname, bool oside) 
{
	afree<FILE*> f(fopen(fname, "rb"), fclose);
	if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
	}


	// map the file indices to array indices
	std::vector<int> facemap;
	std::vector<int> vertmap;


	int linenum=0;
	char delims[] = " \t\n{}=()";
	char line[2048];

	while (fgets(line, 2040, f)) {
		linenum++;

		char *token = strtok(line, delims);

		if (!token || *token=='#') continue;

		if (stricmp(token, "Vertex") == 0)
		{
			token = strtok(NULL, delims);
			if (!token) {
				cerr<<"Mesh::ReadHHM() - error parsing vertex at line "<<linenum<<endl;
				Clear();	return false;
			}

			int vindex = atoi(token);


			value_type ploc[3];

			for (int i=0; i<3; i++) {
				token = strtok(NULL, delims);
				if (!token) {
					cerr<<"Mesh::ReadHHM() - error parsing vertex at line "<<linenum<<endl;
					Clear();	return false;
				}

				ploc[i] = (value_type)atof(token);
			}


			if (oside) {

				while (1) {

					if (!token) {
						cerr<<"Mesh::ReadHHM() - error parsing vertex at line "<<linenum<<endl;
						Clear();	return false;
					}

					if (stricmp(token, "Opos") == 0) {
						for (int i=0; i<3; i++) {
							token = strtok(NULL, delims);
							if (!token) {
								cerr<<"Mesh::ReadHHM() - error parsing vertex at line "<<linenum<<endl;
								Clear();	return false;
							}

							ploc[i] = (value_type)atof(token);
						}

						break;
					}

					token = strtok(NULL, delims);

				}

			}


			if ((int)vertmap.size() <= vindex)
				vertmap.resize(vindex+1);

			vertmap[vindex] = verts.size();
			verts.push_back(TriangleMeshVertex(Point3(ploc[0], ploc[1], ploc[2])));

		} else if (stricmp(token, "Face") == 0) {

			token = strtok(NULL, delims);
			if (!token) {
				cerr<<"Mesh::ReadHHM() - error parsing face at line "<<linenum<<endl;
				Clear();	return false;
			}

			int findex = atoi(token);



			int vi[3];

			for (int i=0; i<3; i++) {
				token = strtok(NULL, delims);
				if (!token) {
					cerr<<"Mesh::ReadHHM() - error parsing face at line "<<linenum<<endl;
					Clear();	return false;
				}

				vi[i] = atoi(token);
			}


			if ((int)facemap.size() <= findex)
				facemap.resize(findex+1);

			facemap[findex] = faces.size();
			faces.push_back(TriangleMeshFace(vi));

		}
	}


	// we've read all the data - build the actual structures now
    build_structures(facemap, vertmap);
	return true;
}

template <class T>
bool tTriangleMesh<T>::read_a_line(FILE* f, char* line, int maxlinelen) 
{
    do
    {
        fgets(line, maxlinelen, f);
    } while (!feof(f) && (line[0] == '#'));
    if (feof(f)) return false;
    else return true;
}


// read a .sma streaming mesh format (to in-core)
template <class T>
bool tTriangleMesh<T>::ReadSMA(const char *fname) 
{
    afree<FILE*> f(fopen(fname, "r"), fclose);
    if (f==0) {
	cerr<<"couldn't open file "<<fname<<endl;
	return false;
    }

    char line[2048];

    while (read_a_line(f, line, 2048)) {
	if (line[0] == 'v') {
	    float x,y,z;
	    sscanf(line, "v %f %f %f", &x, &y, &z);
	    verts.push_back(TriangleMeshVertex(Point3(x, y, z)));
	} else if (line[0] == 'f') {
	    
	    int vi[3];
	    sscanf(line, "f %d %d %d", vi, vi+1, vi+2);
	    vi[0]--;
	    vi[1]--;
	    vi[2]--;
	    faces.push_back(TriangleMeshFace(vi));
	}
    }

    // we've read all the data - build the actual structures now
    std::vector<int> facemap, vertmap;
    IdentityMap(facemap, faces.size());
    IdentityMap(vertmap, verts.size());
    build_structures(facemap, vertmap);

    return true;
}


// read a .off file
template <class T>
bool tTriangleMesh<T>::ReadOFF(const char *fname) 
{
	afree<FILE*> f(fopen(fname, "rb"), fclose);
	if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
	}

	char line[2048];
    /*
     * Read header
     */
    read_a_line(f, line, 2048);
    if (strnicmp(line, "off", 3) != 0)
    {
		cerr<<"Mesh::ReadOFF() - Not a .off file" << endl;
		Clear();	return false;
    }
    read_a_line(f, line, 2048);
    int vnum, fnum;
    sscanf(line, "%d %d", &vnum, &fnum);

    for (int i = 0; i < vnum; ++i)
    {
        if (!read_a_line(f, line, 2048))
        {
		    cerr<<"Mesh::ReadOFF() - EOF" << endl;
		    Clear();	return false;
        }
        float x,y,z;
        sscanf(line, "%f %f %f", &x, &y, &z);
		verts.push_back(TriangleMeshVertex(Point3(x, y, z)));
    }

    for (int i = 0; i < fnum; ++i)
    {
        if (!read_a_line(f, line, 2048))
        {
		    cerr<<"Mesh::ReadOFF() - EOF" << endl;
		    Clear();	return false;
        }
        int n, vi[3];
        sscanf(line, "%d %d %d %d", &n, vi, vi+1, vi+2);
        if (n != 3)
        {
		    cerr<<"Mesh::ReadOFF() - only triangle meshaes are suppored" << endl;
            // continue anyway
        }
		faces.push_back(TriangleMeshFace(vi));
    }

	// we've read all the data - build the actual structures now
	std::vector<int> facemap, vertmap;
	IdentityMap(facemap, fnum);
	IdentityMap(vertmap, vnum);
    build_structures(facemap, vertmap);

	return true;
}

template <class T>
bool tTriangleMesh<T>::WriteHHM(const char *fname) const
{
	afree<FILE*> f(fopen(fname, "wb"), fclose);
	if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
	}

	for (unsigned i=0; i<verts.size(); i++) {
		fprintf(f, "Vertex %d  %g %g %g\n", (i+1), verts[i].point[0], verts[i].point[1], verts[i].point[2]);
	}

	for (unsigned i=0; i<faces.size(); i++) {
		fprintf(f, "Face %d  %d %d %d\n", (i+1), faces[i].verts[0]+1, faces[i].verts[1]+1, faces[i].verts[2]+1);
	}

	return true;
}


template <class T>
bool tTriangleMesh<T>::WriteOFF(const char *fname) const 
{
	afree<FILE*> f(fopen(fname, "wb"), fclose);
	if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
	}

    unsigned NV = verts.size();
    unsigned NF = faces.size();
    fprintf(f,"OFF\n%d %d 0\n", NV, NF);

	for (unsigned i=0; i<NV; i++) {
		fprintf(f, "%g %g %g\n", verts[i].point[0], verts[i].point[1], verts[i].point[2]);
	}

	for (unsigned i=0; i<NF; i++) {
		fprintf(f, "3 %d %d %d\n", faces[i].verts[0], faces[i].verts[1], faces[i].verts[2]);
	}

	return true;
}

template <class T>
bool tTriangleMesh<T>::Write(const char *fname) const
{
    const char* pdot = strrchr(fname, '.');
    if (pdot == 0) return false;
    if (stricmp(pdot+1, "m") == 0) return WriteHHM(fname);
    else if (stricmp(pdot+1, "off") == 0) return WriteOFF(fname);
    else return false;
}

template <class T>
bool tTriangleMesh<T>::Read(const char *fname) 
{
    const char* pdot = strrchr(fname, '.');
    if (pdot == 0) return false;
    if (stricmp(pdot+1, "m") == 0) return ReadHHM(fname);
    else if (stricmp(pdot+1, "off") == 0) return ReadOFF(fname);
    else if (stricmp(pdot+1, "sma") == 0) return ReadSMA(fname);
    else return false;
}

#if 0
DID NOT COMPLETE: using conversion from and to "indexted triangle set"
/*
 * Read a line / skipping comments:
 *  lines that begin with #
 * N - line length
 * Return true if EOF
 */
static bool get_a_line(FILE* f, char* line, int N)
{
    do
    {
        if (fgets(line, N, f) == 0) return true; // EOF
    } while (line[0] == '#');
    return false;
}

template <class T>
bool tTriangleMesh<T>::ReadOFF(const char* fname)
{
	afree<FILE*> f(fopen(fname, "rb"), fclose);
    if (f==0) return false;
	// map the file indices to array indices
	std::map<int, int> facemap;
	std::map<int, int> vertmap;
    char line[1000];
    get_a_line(f, line, 1000);
    if (strnicmp(line, 'OFF', 3) != 0)
    {
        // Not an OFF file
        return false;
    }
    get_a_line(f, line, 1000);
    int V,F;
    sscanf(line, "%d %d", &V, &F);

    build_structures(facemap, vertmap);
}

template <class T>
void tTriangleMesh<T>::WriteOFF(const char* fname)
{
}

#endif // 0

template <class T>
void tTriangleMesh<T>::IdentityMap(std::vector<int> &map, int size) {
	map.resize(size);
	for (int i=0; i<size; i++) {
		map[i] = i;
	}
}


template <class T>
void tTriangleMesh<T>::build_structures(const std::vector<int> &facemap, const std::vector<int> &vertmap)
{
	// convert all the indices to the array indices
	for (typename face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {
		for (int i=0; i<3; i++) {
			f->verts[i] = vertmap[f->verts[i]];
		}
	}

	// set the somefaces
	for (typename vertex_list::iterator v=verts.begin(); v!=verts.end(); ++v) {
		v->someface = -1;
	}
	for (typename face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {
		for (int i=0; i<3; i++) {
			verts[f->verts[i]].someface = FaceIndex(*f);
		}
	}


	// build the adjacency info
	vector< vector<int> > vertfaces(verts.size());
	for (unsigned i=0; i<vertfaces.size(); i++) {
		vertfaces[i].reserve(7);
	}

	for (typename face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {
		for (int i=0; i<3; i++) {
			vertfaces[f->verts[i]].push_back(FaceIndex(*f));
		}
	}

	bool nonmanifold = false;
	for (typename face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {
		for (int i=0; i<3; i++) {

			int v0 = f->verts[i];
			int v1 = f->verts[(i+1)%3];

			// look for a face with the edge v1,v0
			bool found=false;
			for (unsigned vfi=0; vfi<vertfaces[v0].size(); vfi++) {
				int vf = vertfaces[v0][vfi];
				if (faces[vf].EdgeIndexCCW(v1, v0) != -1) {
					f->nbrs[i] = vf;
					if (found) {
						cerr<<"more than one matching triangle found: faces["<<vf<<"]"<<endl;
						nonmanifold = true;
					}
					found=true;
//					break;
				}
			}
		}
	}


	if (nonmanifold) {
	    cerr<<"multi-pass setting normals due to nonmanifoldness"<<endl;

		for (int v=0; v<verts.size(); v++) {
			verts[v].normal = Vector3(0,0,0);
		}

		for (int f=0; f<faces.size(); f++) {
			Vector3 norm = FaceNormal(faces[f]);
			verts[faces[f].verts[0]].normal += norm;
			verts[faces[f].verts[1]].normal += norm;
			verts[faces[f].verts[2]].normal += norm;
		}

		for (int v=0; v<verts.size(); v++) {
			verts[v].normal.normalize();
		}

	} else {
	    SetNormals();
	}


	tModel<T>::invalidate_all();
}


// set the normals
template <class T>
void tTriangleMesh<T>::SetNormals() {
	for (typename vertex_list::iterator v=verts.begin(); v!=verts.end(); ++v) {
		v->normal = Vector3(0,0,0);
		for (VertexFaceIterator f(*this, VertexIndex(*v)); !f.done(); ++f) {

			Vector3 e1 = verts[(*f).verts[1]].point - verts[(*f).verts[0]].point;
			Vector3 e2 = verts[(*f).verts[2]].point - verts[(*f).verts[0]].point;
			Vector3 fn = e1.cross(e2);
			if (fn.length() == 0) {
				cerr<<"skipping normal of degenerate face "<<FaceIndex(*f)<<endl;
			} else {
				fn.normalize();
				v->normal += FaceNormal(*f);
			}
		}
		v->normal.normalize();
	}
}



template <class T>
void tTriangleMesh<T>::Clear() {
	verts.clear();
	faces.clear();
}



// Model stuff
template <class T>
void tTriangleMesh<T>::compute_bounding_box() const {

        tModel<T>::_bounding_box = Box3(verts[0].point, verts[0].point);
	for (typename vertex_list::const_iterator v=verts.begin(); v!=verts.end(); ++v) {
	    tModel<T>::_bounding_box.update(v->point);
	}
	tModel<T>::_is_bounding_box_valid = true;
}


template <class T>
void tTriangleMesh<T>::compute_centroid() const {

	tModel<T>::_centroid = Point3(0,0,0);

	for (typename vertex_list::const_iterator v=verts.begin(); v!=verts.end(); ++v) {
		tModel<T>::_centroid.add(v->point);
	}
	tModel<T>::_centroid.scalar_scale(1.0 / verts.size());
	tModel<T>::_is_centroid_valid = true;
}




// iterators

template<class T>
tTriangleMesh<T>::VertexFaceIterator::VertexFaceIterator(tTriangleMesh<T> &_mesh, int vi) {

	first=true;
	mesh = &_mesh;
	vertex = vi;

	// rotate clockwise as far as possible
	cur_face = mesh->verts[vi].someface;
	if (cur_face != -1) {
		do {
			TriangleMeshFace &f = mesh->faces[cur_face];
			int cindex = f.VertIndex(vi);
			if (f.nbrs[cindex] == -1) break;
			cur_face = f.nbrs[cindex];
		} while (cur_face != mesh->verts[vi].someface);
		first_face = cur_face;
	}
}

template<class T>
bool tTriangleMesh<T>::VertexFaceIterator::done() {
	if (cur_face == -1) return true;
	return false;
}


template<class T>
TriangleMeshFace& tTriangleMesh<T>::VertexFaceIterator::operator*() {
	return mesh->faces[cur_face];
}


template<class T>
typename tTriangleMesh<T>::VertexFaceIterator& tTriangleMesh<T>::VertexFaceIterator::operator++() {

	TriangleMeshFace &f = mesh->faces[cur_face];
	int nindex = (f.VertIndex(vertex) + 2) % 3;
	cur_face = f.nbrs[nindex];

	if (cur_face==first_face)
		cur_face = -1;

	first = false;

	return *this;
}


template<class T>
tTriangleMesh<T>::VertexFaceIteratorI::VertexFaceIteratorI(const tTriangleMesh<T> &_mesh, int vi) {

	first=true;
	mesh = &_mesh;
	vertex = vi;

	// rotate clockwise as far as possible
	cur_face = mesh->verts[vi].someface;
	if (cur_face != -1) {
		do {
			const TriangleMeshFace &f = mesh->faces[cur_face];
			int cindex = f.VertIndex(vi);
			if (f.nbrs[cindex] == -1) break;
			cur_face = f.nbrs[cindex];
		} while (cur_face != mesh->verts[vi].someface);
		first_face = cur_face;
	}
}

template<class T>
bool tTriangleMesh<T>::VertexFaceIteratorI::done() {
	if (cur_face == -1) return true;
	return false;
}


template<class T>
int tTriangleMesh<T>::VertexFaceIteratorI::operator*() {
	return cur_face;
}


template<class T>
typename tTriangleMesh<T>::VertexFaceIteratorI& tTriangleMesh<T>::VertexFaceIteratorI::operator++() {

	const TriangleMeshFace &f = mesh->faces[cur_face];
	int nindex = (f.VertIndex(vertex) + 2) % 3;
	cur_face = f.nbrs[nindex];

	if (cur_face==first_face)
		cur_face = -1;

	first = false;

	return *this;
}





template<class T>
tTriangleMesh<T>::VertexVertexIterator::VertexVertexIterator(tTriangleMesh<T> &_mesh, int vi) : vfi(_mesh, vi) {
	first = true;
	int v = (*vfi).VertIndex(vi);
	onedge = ((*vfi).nbrs[v] == -1);
}


template<class T>
bool tTriangleMesh<T>::VertexVertexIterator::done() {
	return vfi.done();
}


template<class T>
typename tTriangleMesh<T>::TriangleMeshVertex& tTriangleMesh<T>::VertexVertexIterator::operator*() {

	TriangleMeshFace &f = *vfi;
	int vi = f.VertIndex(vfi.vertex);

	if (onedge && first) {
		return vfi.mesh->verts[f.verts[(vi+1)%3]];
	} else {
		return vfi.mesh->verts[f.verts[(vi+2)%3]];
	}
}


template<class T>
typename tTriangleMesh<T>::VertexVertexIterator& tTriangleMesh<T>::VertexVertexIterator::operator++() {

	if (!(onedge && first))
		++vfi;
	first = false;

	return *this;
}




template<class T>
tTriangleMesh<T>::VertexVertexIteratorI::VertexVertexIteratorI(const tTriangleMesh<T> &_mesh, int vi) : vfi(_mesh, vi) {
	first = true;
	int ff = *vfi;
	int v = _mesh.faces[ff].VertIndex(vi);
	onedge = (_mesh.faces[ff].nbrs[v] == -1);
}


template<class T>
bool tTriangleMesh<T>::VertexVertexIteratorI::done() {
	return vfi.done();
}


template<class T>
int tTriangleMesh<T>::VertexVertexIteratorI::operator*() {

	const TriangleMeshFace &f = vfi.mesh->faces[*vfi];
	int vi = f.VertIndex(vfi.vertex);

	if (onedge && first) {
		return f.verts[(vi+1)%3];
	} else {
		return f.verts[(vi+2)%3];
	}
}


template<class T>
typename tTriangleMesh<T>::VertexVertexIteratorI& tTriangleMesh<T>::VertexVertexIteratorI::operator++() {

	if (!(onedge && first))
		++vfi;
	first = false;

	return *this;
}




// do one iteration of loop subdivision
template<class T>
void tTriangleMesh<T>::LoopSubdivide(bool keep_current_positions) {

	// we'll put all the verts/faces into an its, then convert back to setup the adjacency info
	tIndexedTriangleSet<T> its;


	// find the new positions and update all the existing points
	for (typename vertex_list::iterator v=verts.begin(); v!=verts.end(); ++v) {

		vector<int> ring;

		VertexVertexIteratorI vi(*this,  VertexIndex(*v));
		for ( ; !vi.done(); ++vi) {
			ring.push_back(*vi);
		}

		Point3 np(0,0,0);

		if (keep_current_positions) {
		  np = v->point;
		} else {
		  if (vi.isonedge()) {
		    np.add_scaled(v->point, 6.0/8.0);
		    np.add_scaled(verts[ring.front()].point, 1.0/8.0);
		    np.add_scaled(verts[ring.back()].point, 1.0/8.0);
		  } else {
		    //			value_type alpha = (3 + 2*cos(2*M_PI / ring.size()));
		    //			alpha = (40 - alpha*alpha) / 64;
		    //			alpha /= ring.size();
		    double alpha = (ring.size()>3) ? (3.0/(8.0*ring.size())) : (3.0/16.0);	// warren weights
		    
		    np.add_scaled(v->point, 1.0-ring.size()*alpha);
		    
		    for (unsigned i=0; i<ring.size(); i++) {
		      np.add_scaled(verts[ring[i]].point, alpha);
		    }
		  }
		}

		its.add_vertex(np);
	}


	// create and set the new vertex positions
	int nnverts = 0;
	vector<int> nvi(faces.size()*3);

	for (face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {

		int fi = FaceIndex(*f);

		for (int i=0; i<3; i++) {
			if (f->nbrs[i] < fi) {

				Point3 np(0,0,0);
				if (keep_current_positions || f->nbrs[i] < 0) {
					np.add_scaled(verts[f->verts[i]].point, 0.5);
					np.add_scaled(verts[f->verts[(i+1)%3]].point, 0.5);
				} else {
					np.add_scaled(verts[f->verts[i]].point, 3.0/8.0);
					np.add_scaled(verts[f->verts[(i+1)%3]].point, 3.0/8.0);

					const TriangleMeshFace &nf = faces[f->nbrs[i]];
					np.add_scaled(verts[f->verts[(i+2)%3]].point, 1.0/8.0);
					np.add_scaled(verts[nf.verts[(nf.VertIndex(f->verts[i])+1)%3]].point, 1.0/8.0);
				}

				its.add_vertex(np);
				nvi[fi*3 + i] = verts.size() + nnverts++;
			}
			else
				nvi[fi*3 + i] = -1;
		}
	}

	for (int fi=0; fi<(int)faces.size(); fi++) {

		const TriangleMeshFace &mf = faces[fi];

		int aoeu[3] = { nvi[fi*3 + 0], nvi[fi*3 + 1], nvi[fi*3 + 2] };

		for (int i=0; i<3; i++) {
			int aoeu2 = faces[fi].nbrs[i];
			if (faces[fi].nbrs[i] > fi) {
				if (nvi[fi*3 + i] >= 0)
					BREAK;
			} else {
				if (nvi[fi*3 + i] < 0)
					BREAK;
			}
		}

	}



	// setup all the new faces
	for (face_list::iterator f=faces.begin(); f!=faces.end(); ++f) {

		int fi = FaceIndex(*f);

		int tnvi[3];
		for (int i=0; i<3; i++) {
			if (f->nbrs[i] < fi) {
				tnvi[i] = nvi[fi*3 + i];
			} else {

				const TriangleMeshFace &nf = faces[f->nbrs[i]];

				int aoeu[3] = { nvi[f->nbrs[i]*3 + 0], nvi[f->nbrs[i]*3 + 1], nvi[f->nbrs[i]*3 + 2] };

				int ni = (nf.VertIndex(f->verts[i]) + 2) % 3;
				tnvi[i] = nvi[f->nbrs[i]*3 + ni];
			}
		}

		its.insert_triangle((unsigned)f->verts[0], (unsigned)tnvi[0], (unsigned)tnvi[2], false, 0.0);
		its.insert_triangle((unsigned)f->verts[1], (unsigned)tnvi[1], (unsigned)tnvi[0], false, 0.0);
		its.insert_triangle((unsigned)f->verts[2], (unsigned)tnvi[2], (unsigned)tnvi[1], false, 0.0);
		its.insert_triangle((unsigned)tnvi[0], (unsigned)tnvi[1], (unsigned)tnvi[2], false, 0.0);
	}

	Clear();

	ITS2TM(its, *this);
}




// create the triangle strip structure for rendering
// adapted from:
/*
Szymon Rusinkiewicz
Princeton University

TriMesh_tstrips.cc
Code for dealing with triangle strips
*/
template<class T>
void tTriangleMesh<T>::BuildStrips(int maxlen) {

	tstrips.clear();

	std::vector<int> todo;
	std::vector<signed char> face_avail(faces.size());
	unsigned i;
	for (i = 0; i < faces.size(); i++) {
		face_avail[i] = (faces[i].nbrs[0] != -1) +
						(faces[i].nbrs[1] != -1) +
						(faces[i].nbrs[2] != -1);
		if (face_avail[i] == 1)
			todo.push_back(i);
	}

	tstrips.reserve(faces.size() * 2);

	int nstrips = 0;
	i = 0;
	while (i < faces.size()) {
		int next;
		if (todo.empty()) {
			next = i++;
		} else {
			next = todo.back();
			todo.pop_back();
		}
		if (face_avail[next] < 0)
			continue;
		BuildStrip(next, face_avail, todo, maxlen);
		nstrips++;
	}


	// Collect all single triangles to be at the end of the list of tstrips
	{
		vector<int> tris;
		int n = 0, offset = 0;
		bool have_tri = false, bad_strip = false;
		for (int i = 0; i < (int)tstrips.size(); i++) {
			if (n == 0) {
				n = tstrips[i];
				bad_strip = (n < 3);
				have_tri = (n == 3);
				n++;
			}
			if (bad_strip) {
				offset++;
			} else if (have_tri) {
				tris.push_back(tstrips[i]);
				offset++;
			} else if (offset > 0) {
				tstrips[i-offset] = tstrips[i];
			}
			n--;
		}
		if (offset == 0)
			return;

		tstrips.erase(tstrips.end() - offset, tstrips.end());
		tstrips.insert(tstrips.end(), tris.begin(), tris.end());
	}

	// renumber the vertices so they're more cache-coherent
	if (0){
		vector<int> vertorder(tstrips.size());
		vertorder.clear();

		{
			const int *index = &tstrips[0];
			const int *end = &(*tstrips.end());
			while (index<end) {
				int len = *index;
				index++;
				for (int i=0; i<len; i++) {
					vertorder.push_back(index[i]);
				}
				index += len;
			}
		}

		vector<int> old_to_new(verts.size());
		vector<int> new_to_old(verts.size());
		new_to_old.clear();


		for (unsigned i=0; i<old_to_new.size(); i++) {
			old_to_new[i] = -1;
		}


		for (unsigned i=0; i<vertorder.size(); i++) {
			if (old_to_new[vertorder[i]] < 0) {
				old_to_new[vertorder[i]] = new_to_old.size();
				new_to_old.push_back(vertorder[i]);
			}
		}

		vertex_list	nverts(verts.size());
		for (unsigned i=0; i<verts.size(); i++) {
			nverts[old_to_new[i]] = verts[i];
		}
		verts = nverts;
    
		for (unsigned i=0; i<faces.size(); i++) {
			for (int j=0; j<3; j++) {
				faces[i].verts[j] = old_to_new[faces[i].verts[j]];
			}
		}


		{
			int *index = &tstrips[0];
			int *end = &(*tstrips.end());
			while (index<end) {
				int len = *index;
				index++;
				for (int i=0; i<len; i++) {
					index[i] = old_to_new[index[i]];
				}
				index += len;
			}
		}

	}



}

template<class T>
void tTriangleMesh<T>::BuildStrip(int f, std::vector<signed char> &face_avail, std::vector<int> &todo, int maxlen) {
	TriangleMeshFace &face = faces[f];
	if (face_avail[f] == 0) {
		tstrips.push_back(3);
		tstrips.push_back(face.verts[0]);
		tstrips.push_back(face.verts[1]);
		tstrips.push_back(face.verts[2]);
		face_avail[f] = -1;
		return;
	}

	int score[3];
	for (int i = 0; i < 3; i++) {
		score[i] = 0;
		int ae = face.nbrs[(i+1)%3];
		if (ae == -1 || face_avail[ae] < 0)
			continue;
		score[i]++;
		int next_edge = (faces[ae].VertIndex(face.verts[(i+1)%3])+1)%3;
		int nae = faces[ae].nbrs[next_edge];
		if (nae == -1 || face_avail[nae] < 0)
			continue;
		score[i]++;
		if (face_avail[ae] == 2)
			score[i]++;
	}

	int best_score = max3(score[0], score[1], score[2]);
	int best = (score[0] == best_score) ? 0 :
		   (score[1] == best_score) ? 1 : 2;

    int vlast2 = face.verts[ best     ];
	int vlast1 = face.verts[(best+1)%3];
	int vnext  = face.verts[(best+2)%3];
	int dir = 1;

	unsigned tsi=tstrips.size();
	tstrips.push_back(-1);
	tstrips.push_back(vlast2);
	tstrips.push_back(vlast1);

	while (1) {
		tstrips.push_back(vnext);
		face_avail[f] = -1;
		for (int j = 0; j < 3; j++) {
			int ae = faces[f].nbrs[j];
			if (ae == -1)
				continue;
			if (face_avail[ae] > 0)
				face_avail[ae]--;
			if (face_avail[ae] == 1)
				todo.push_back(ae);
		}

		f = faces[f].nbrs[(faces[f].VertIndex(vlast2)+1)%3];
		if (f == -1 || face_avail[f] < 0)
			break;

		if (tstrips.size()-tsi > (unsigned)maxlen)
			break;

		vlast2 = vlast1;
		vlast1 = vnext;
		vnext = faces[f].verts[(faces[f].VertIndex(vlast2)+3+dir)%3];
		dir = -dir;
	}

	tstrips[tsi] = tstrips.size() - tsi - 1;
}


// force instantiation
template class tTriangleMesh<float>;
template class tTriangleMesh<double>;

/*------------------------ Tools --------------------------*/

/*
 * Triangle mesh ---> Indext Triangle Set
 */
template <class T>
void TM2ITS(const tTriangleMesh<T>& tm, tIndexedTriangleSet<T>& its)
{
    its.clear();
    unsigned V = tm.num_vertices();
    for (unsigned i = 0; i < V; ++i)
    {
        its.insert_vertex(tm.vertex(i).point);
    }

    unsigned F = tm.num_faces();
    for (unsigned i = 0; i < F; ++i)
    {
        const TriangleMeshFace& triangle = tm.face(i);
        its.insert_triangle((unsigned)triangle.verts[0], (unsigned)triangle.verts[1], (unsigned)triangle.verts[2], false, 0.0);
    }
}

// Instantiate
template void TM2ITS(const tTriangleMesh<float>& tm, tIndexedTriangleSet<float>& its);
template void TM2ITS(const tTriangleMesh<double>& tm, tIndexedTriangleSet<double>& its);

template <class T>
void ITS2TM(const tIndexedTriangleSet<T>& its, tTriangleMesh<T>& tm)
{
    tm.Clear();

    unsigned V = its.num_vertices();
    for (unsigned i = 0; i < V; ++i)
    {
        tm.verts.push_back(typename tTriangleMesh<T>::TriangleMeshVertex(its.vertex(i)));
    }

    unsigned F = its.num_triangles();
    for (unsigned i = 0; i < F; ++i)
    {
        TriangleMeshFace fi(its.indexed_triangle(i).get_indices());
        tm.faces.push_back(fi);
    }

	std::vector<int> facemap, vertmap;
	tm.IdentityMap(facemap, F);
	tm.IdentityMap(vertmap, V);
    tm.build_structures(facemap, vertmap);
}

// Instantiate
template void ITS2TM(const tIndexedTriangleSet<float>& its, tTriangleMesh<float>& tm);
template void ITS2TM(const tIndexedTriangleSet<double>& its, tTriangleMesh<double>& tm);

GTB_END_NAMESPACE

