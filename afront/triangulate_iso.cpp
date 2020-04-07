
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
#include "triangulate_iso.h"
#include "triangulate_mesh.h"
#include "triangulate_tet.h"

#include "interval.h"


#include <lib/mlslib/NR/nr.h>	// for zbrent?!?!?


const real_type theta_step = M_PI_2/100;

real_type iso_scale = 1;

bool reuse_guidance = false;
bool reuse_seeds = false;

bool endswith(const char *s, const char *end);


#define ISO_CIRCULAR_PROJECTION

extern int curvature_sub;


#define DIM 256

//#define ZORDER_DATA


template <typename T>
int round_to_minus_inf(const T &x) {
    return ((x>=0) ?
			(int)x :
			(int)x - 1);
}

template <typename T>
RegularVolume<T>::RegularVolume() : kdGetPoint(gfpoints)
{
    data = NULL;
    kdtree = NULL;
}


template <typename T>
RegularVolume<T>::~RegularVolume() {
    if (data) delete [] data;
    if (kdtree) delete kdtree;
}


template <typename T>
real_type RegularVolume<T>::GetValue(int x, int y, int z) const {
#ifdef ZORDER_DATA
	int64 i;
	ZCurve::xyz2i(x,y,z,i);
	return (real_type)data[i] * iso_scale;
#else
    return (real_type)data[(z*dim[0]*dim[1] + y*dim[0] + x)] * iso_scale;
#endif
}


template <typename T>
void RegularVolume<T>::Gather(const int cell[3], double nbrs[4][4][4]) const {

    for (int z=0; z<4; z++) {
		int zi = clamp(cell[2]+z-1, 0, dim[2]-1);

		for (int y=0; y<4; y++) {
			int yi = clamp(cell[1]+y-1, 0, dim[1]-1);

			for (int x=0; x<4; x++) {
				int xi = clamp(cell[0]+x-1, 0, dim[0]-1);

#ifdef ZORDER_DATA
				int64 i;
				ZCurve::xyz2i(xi,yi,zi,i);
				nbrs[x][y][z] = (double)data[i] * iso_scale;
#else
				nbrs[x][y][z] = (double)data[(zi*dim[0]*dim[1] + yi*dim[0] + xi)] * iso_scale;
#endif
			}
		}
    }
}


template <typename T>
int RegularVolume<T>::NumBlocks() const {
    return 1;
}

template <typename T>
void RegularVolume<T>::SetActiveBlock(int b) const {
}

template <typename T>
int RegularVolume<T>::GetActiveBlock() const {
    return 0;
}

template <typename T>
int* RegularVolume<T>::GetActiveBlock3(int i3[3]) const {
    i3[0] = 0;
    i3[1] = 0;
    i3[2] = 0;
	return i3;
}

template <typename T>
void RegularVolume<T>::GetBlockBox(int min[3], int max[3]) const {
    for (int i=0; i<3; i++) {
		min[i] = 0;
		max[i] = dim[i];
    }
}

template <typename T>
void RegularVolume<T>::GetBlockBox(int block[3], int min[3], int max[3]) const {
    return GetBlockBox(min, max);
}

template <typename T>
void RegularVolume<T>::GetBlockMinMax(real_type &min, real_type &max) const {
    min = (real_type)(this->min) * iso_scale;
    max = (real_type)(this->max) * iso_scale;
}



template <typename T>
void RegularVolume<T>::WriteValue(int x, int y, int z, FILE *f) const {
#ifdef ZORDER_DATA
	int64 i;
	ZCurve::xyz2i(x,y,z,i);
    fwrite(&data[i], 1, sizeof(T), f);
#else
    fwrite(&data[(z*dim[0]*dim[1] + y*dim[0] + x)], 1, sizeof(T), f);
#endif
}


#ifdef HAS_ZLIB
template <typename T>
void RegularVolume<T>::WriteValue(int x, int y, int z, gzFile f) const {
#ifdef ZORDER_DATA
	int64 i;
	ZCurve::xyz2i(x,y,z,i);
    gzwrite(f, &data[i], sizeof(T));
#else
    gzwrite(f, &data[(z*dim[0]*dim[1] + y*dim[0] + x)], sizeof(T));
#endif
}
#endif

template <>
const char* RegularVolume<unsigned char>::TypeString() const {
    return "uchar";
}
template <>
const char* RegularVolume<char>::TypeString() const {
    return "signed char";
}
template <>
const char* RegularVolume<float>::TypeString() const {
    return "float";
}
template <>
const char* RegularVolume<short>::TypeString() const {
    return "short";
}


template <typename T>
void RegularVolume<T>::SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal) {

    gfpoints = pts;
    ideal_lengths = ideal;

	// cant make empty trees!
	if (pts.size() > 0) {
		kdtree = new kdtree_type(10, Box3::bounding_box(pts), kdGetPoint);
		for (unsigned i=0; i<gfpoints.size(); i++)
			kdtree->Insert(i);
		kdtree->MakeTree();
	} else {
		kdtree = NULL;
	}
}

template <typename T>
void RegularVolume<T>::SetGuidanceFieldDone() {
}


template <typename T>
void RegularVolume<T>::SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds) {

    boundary_pts = bnd_pts;
    boundary_norms = bnd_norms;
    boundary_loops = loops;
    seed_pts = cc_seeds;
}

template <typename T>
void RegularVolume<T>::GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const {

    bnd_pts = boundary_pts;
    bnd_norms = boundary_norms;
    loops = boundary_loops;
    cc_seeds = seed_pts;
}


template <typename T>
void RegularVolume<T>::GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {

    int bmins[3], bmaxs[3];
    GetBlockBox(bmins, bmaxs);

    box = Box3(Point3(1e34,1e34,1e34), Point3(-1e34,-1e34,-1e34));
    for (int i=0; i<3; i++) {
		if (GetAspect(i) > 0) {
			box.set_min(i, bmins[i]*GetAspect(i) + translation[i]);
			box.set_max(i, bmaxs[i]*GetAspect(i) + translation[i]);
		} else {
			box.set_max(i, bmins[i]*GetAspect(i) + translation[i]);
			box.set_min(i, bmaxs[i]*GetAspect(i) + translation[i]);
		}
    }

    pts = boundary_pts;
    norms = boundary_norms;
    loops = boundary_loops;
    seeds = seed_pts;
}


template <typename T>
bool RegularVolume<T>::GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
    return false;
}

template <typename T>
int RegularVolume<T>::GetBlockForPoint(const Point3 &p) const {
	return 0;
}

template <typename T>
void* RegularVolume<T>::OrderedPointTraverseStart(const Point3 &p) {

    if (kdtree) {
		return new typename kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
	} else {
		return NULL;
	}
}


template <typename T>
bool RegularVolume<T>::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {

    if (!ctx) {
		return false;
	}

	typename kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (typename kdtree_type::OrderedIncrementalTraverse*)ctx;
	if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = gfpoints[n];
    ideal = ideal_lengths[n];
    return true;
}


template <typename T>
void RegularVolume<T>::OrderedPointTraverseEnd(void *ctx) {
	if (ctx) {
		typename kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (typename kdtree_type::OrderedIncrementalTraverse*)ctx;
		delete kdOrderedTraverse;
	}
}


template <typename T>
bool RegularVolume<T>::ReadSource(const char *fname, bool bigendian) {

    if (endswith(fname, ".gz")) {

#ifdef HAS_ZLIB
		gtb::afree<gzFile> f(gzopen(fname, "rb"), gzclose);

		if (f==0) {
			cerr<<"couldn't open file "<<fname<<endl;
			return false;
		}


		if (dim[2] == 0) {
			cerr<<"dim[2] not set for gz file!"<<endl;
			return false;
		}

#ifdef ZORDER_DATA
		int zdim = std::max(std::max(dim[0],dim[1]), dim[2]);
		int i=1<<30;
		while ((zdim & i) == 0) {
			i>>=1;
		}

		if (zdim != i)
			zdim = i<<1;

		data = new T[zdim*zdim*zdim];
#else
		data = new T[dim[0]*dim[1]*dim[2]];
#endif

		bool first=true;
		for (int z=0; z<dim[2]; z++) {
			for (int y=0; y<dim[1]; y++) {
				for (int x=0; x<dim[0]; x++) {

					T s;
					gzread(f, &s, sizeof(T));

					if (bigendian) {
						char* endianswap = (char*)&s;
						for (int j=0; j<sizeof(T)/2; j++) {
							std::swap(endianswap[j], endianswap[sizeof(T)-j-1]);
						}
					}

#ifdef ZORDER_DATA
					int64 i;
					ZCurve::xyz2i(x,y,z,i);
					data[i] = s;
#else
					data[(z*dim[0]*dim[1] + y*dim[0] + x)] = s;
#endif
					if (first) {
						first = false;
						min = s;
						max = s;
					} else {
						min = std::min(min, s);
						max = std::max(max, s);
					}
				}
			}
		}

#else

		cerr<<"gzip files not supported!"<<endl;
		return false;

#endif

    } else {

		gtb::afree<FILE*> f(fopen(fname, "rb"), fclose);

		if (f==0) {
			cerr<<"couldn't open file "<<fname<<endl;
			return false;
		}


		if (dim[2] == 0) {

			fseek(f, 0, SEEK_END);
			int size = ftell(f);

			if (size % (sizeof(T)*dim[0]*dim[1])) {
				cerr<<"volume file size strange: "<<size<<"!"<<endl;
				return false;
			}
	  
			dim[2] = size/(sizeof(T)*dim[0]*dim[1]);
		}


#ifdef ZORDER_DATA
		int zdim = std::max(std::max(dim[0],dim[1]), dim[2]);
		int i=1<<30;
		while ((zdim & i) == 0) {
			i>>=1;
		}

		if (zdim != i)
			zdim = i<<1;

		data = new T[zdim*zdim*zdim];
#else
		data = new T[dim[0]*dim[1]*dim[2]];
#endif


		bool first=true;

		fseek(f, 0, SEEK_SET);

		for (int z=0; z<dim[2]; z++) {
			for (int y=0; y<dim[1]; y++) {
				for (int x=0; x<dim[0]; x++) {

					T s;
					fread(&s, sizeof(T), 1, f);

					if (bigendian) {
						char* endianswap = (char*)&s;
						for (int j=0; j<sizeof(T)/2; j++) {
							std::swap(endianswap[j], endianswap[sizeof(T)-j-1]);
						}
					}

#ifdef ZORDER_DATA
					int64 i;
					ZCurve::xyz2i(x,y,z,i);
					data[i] = s;
#else
					data[(z*dim[0]*dim[1] + y*dim[0] + x)] = s;
#endif

					if (first) {
						first = false;
						min = s;
						max = s;
					} else {
						min = std::min(min, s);
						max = std::max(max, s);
					}
				}
			}
		}
    }

    return true;
}


template <typename T>
bool RegularVolume<T>::ReadNRRD(const NHDR &nhdr) {

    aspect[0] = nhdr.spacings[0];
    aspect[1] = nhdr.spacings[1];
    aspect[2] = nhdr.spacings[2];

    dim[0] = nhdr.sizes[0];
    dim[1] = nhdr.sizes[1];
    dim[2] = nhdr.sizes[2];

    translation = nhdr.translation;

    return ReadSource(nhdr.datafile, nhdr.endian);
}


template <typename T>
void RegularVolume<T>::FromSegmentation(const Volume &v, int mat) {

	aspect[0] = v.GetAspect(0);
	aspect[1] = v.GetAspect(1);
	aspect[2] = v.GetAspect(2);

	dim[0] = v.GetDim(0);
	dim[1] = v.GetDim(1);
	dim[2] = v.GetDim(2);

	translation = v.GetTranslation();

	spline = v.spline;

	data = new T[dim[0]*dim[1]*dim[2]];
	min = -1;
	max = 1;

	for (int z=0; z<dim[2]; z++) {
		for (int y=0; y<dim[1]; y++) {
			for (int x=0; x<dim[0]; x++) {
				int i = (z*dim[0]*dim[1] + y*dim[0] + x);

				data[i] = (T)0;
				real_type vv = v.GetValue(x,y,z);
				if (vv == mat) {
					data[i] = (T)1;
				}
			}
		}
	}
}


template <typename T>
void RegularVolume<T>::FromMultiMaterial(vector<Volume*> &vols, int mat) {

	aspect[0] = vols[0]->GetAspect(0);
	aspect[1] = vols[0]->GetAspect(1);
	aspect[2] = vols[0]->GetAspect(2);

	dim[0] = vols[0]->GetDim(0);
	dim[1] = vols[0]->GetDim(1);
	dim[2] = vols[0]->GetDim(2);

	translation = vols[0]->GetTranslation();

	spline = vols[0]->spline;

	data = new T[dim[0]*dim[1]*dim[2]];
	min = (T)INFINITY;
	max = (T)-INFINITY;

	for (int z=0; z<dim[2]; z++) {
		for (int y=0; y<dim[1]; y++) {
			for (int x=0; x<dim[0]; x++) {
				
				real_type vmax = -INFINITY;
				real_type vmin = INFINITY;
				for (int i=0; i<vols.size(); i++) { 
					if (i==mat) continue;
					vmax = std::max(vmax, vols[i]->GetValue(x,y,z));
					vmin = std::min(vmin, vols[i]->GetValue(x,y,z));
				}

				data[(z*dim[0]*dim[1] + y*dim[0] + x)] = vols[mat]->GetValue(x,y,z)-vmax;

				min = std::min(min, data[(z*dim[0]*dim[1] + y*dim[0] + x)]);
				max = std::max(max, data[(z*dim[0]*dim[1] + y*dim[0] + x)]);
			}
		}
	}
}



template <typename T>
void RegularVolume<T>::Smoothed(const Volume &v, int width) {

	aspect[0] = v.GetAspect(0);
	aspect[1] = v.GetAspect(1);
	aspect[2] = v.GetAspect(2);

	dim[0] = v.GetDim(0);
	dim[1] = v.GetDim(1);
	dim[2] = v.GetDim(2);

	translation = v.GetTranslation();

	spline = v.spline;

	data = new T[dim[0]*dim[1]*dim[2]];

	for (int z=0; z<dim[2]; z++) {
		for (int y=0; y<dim[1]; y++) {
			for (int x=0; x<dim[0]; x++) {
				int i = (z*dim[0]*dim[1] + y*dim[0] + x);

				data[i] = (T)0;
				real_type w = 0;

				for (int z2=z-width; z2<=z+width; z2++) {
					for (int y2=y-width; y2<=y+width; y2++) {
						for (int x2=x-width; x2<=x+width; x2++) {
							data[i] += v.GetValue(clamp(x2,0,dim[0]-1),
												  clamp(y2,0,dim[1]-1),
												  clamp(z2,0,dim[2]-1));
 							w += 1;
						}
					}
				}

				data[i] /= w;

				min = std::min(min, data[i]);
				max = std::max(max, data[i]);
			}
		}
	}

	min = data[0];
	max = data[0];

	for (int i=1; i<dim[0]*dim[1]*dim[2]; i++) {
		min = std::min(min, data[i]);
		max = std::max(max, data[i]);
	}
}

template class RegularVolume<char>;
template class RegularVolume<float>;



bool ReadNHDR(const char *fname, NHDR &nhdr) {
    gtb::afree<FILE*> f(fopen(fname, "rb"), fclose);

    if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
    }

    char sst[1024];
    nhdr.type = NHDR::Unknown;

    nhdr.sizes[0] = -1;
    nhdr.sizes[1] = -1;
    nhdr.sizes[2] = -1;

    float faspect[3];
    nhdr.spacings[0] = 1;
    nhdr.spacings[1] = 1;
    nhdr.spacings[2] = 1;

    char datafile[1024];
    datafile[0]=0;
    nhdr.datafile[0]=0;

    char endians[1024];
    nhdr.endian = false;

    float forigin[3] = {0,0,0};
    nhdr.translation[0] = 0;
    nhdr.translation[1] = 0;
    nhdr.translation[2] = 0;


    char line[1024];
    while (fgets(line, 1000, f)) {
		if (1 == sscanf(line, "type: %s", sst)) {
			if (!stricmp(sst, "float")) nhdr.type=NHDR::Float;
			if (!stricmp(sst, "uchar")) nhdr.type=NHDR::UChar;
			if (!stricmp(line+6, "unsigned")) nhdr.type=NHDR::UChar;
			if (!stricmp(line+6, "unsigned char\n")) nhdr.type=NHDR::UChar;
			if (!stricmp(line+6, "unsigned short\n")) nhdr.type=NHDR::Short;
			if (!stricmp(sst, "short")) nhdr.type=NHDR::Short;
		} else if (3 == sscanf(line, "sizes: %d %d %d", &nhdr.sizes[0], &nhdr.sizes[1], &nhdr.sizes[2])) {
		} else if (3 == sscanf(line, "spacings: %g %g %g", &faspect[0], &faspect[1], &faspect[2])) {
			nhdr.spacings[0]=faspect[0]; nhdr.spacings[1]=faspect[1]; nhdr.spacings[2]=faspect[2]; 
		} else if (1 == sscanf(line, "data file: %s", datafile)) {
		} else if (1 == sscanf(line, "endian: %s", endians)) {
			nhdr.endian = (!stricmp(endians, "big"));
		} else if (3 == sscanf(line, "space origin: (%g,%g,%g)", &forigin[0], &forigin[1], &forigin[2])) {
			nhdr.translation = Vector3(forigin[0], forigin[1], forigin[2]);
		}
    }

    if (nhdr.type == NHDR::Unknown ||
		nhdr.sizes[0]<0 || nhdr.sizes[1]<0 || nhdr.sizes[2]<0 ||
		nhdr.spacings[0]==0 || nhdr.spacings[1]==0 || nhdr.spacings[2]==0 ||
		datafile[0]==0) {
		cerr<<"error reading nrrd header"<<endl;
		cerr<<"type: "<<nhdr.type<<endl;
		cerr<<"sizes: "<<nhdr.sizes[0]<<" "<<nhdr.sizes[1]<<" "<<nhdr.sizes[2]<<endl;
		cerr<<"spacings: "<<nhdr.spacings[0]<<" "<<nhdr.spacings[1]<<" "<<nhdr.spacings[2]<<endl;
		cerr<<"datafile: "<<datafile<<endl;
		return false;
    }

    strcpy(nhdr.datafile, fname);
    char *slash = strrchr(nhdr.datafile, '/');

    if (!slash) {
		strcpy(nhdr.datafile, datafile);
    } else {
		strcpy(slash+1, datafile);
    }

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////


OOCVolume::OOCVolume() {
    blocksize[0] = 0;
    blocksize[1] = 0;
    blocksize[2] = 0;
    nblocks[0] = 0;
    nblocks[1] = 0;
    nblocks[2] = 0;
    active_index = -1;
    active_block[0] = -1;
    active_block[1] = -1;
    active_block[2] = -1;
    active_zindex = -1;
    blockToLoadASAP = -1;
	num_block_loads = 0;
    active_region_loading = false;
    guidance_computed = reuse_guidance;
    load_counter = 0;

    prefetcher = NULL;
#ifndef WIN32
	pthread_mutexattr_t attr;
	pthread_mutexattr_init(&attr);
	pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE_NP);
	pthread_mutex_init(&fetch_mutex, &attr);

	pthread_cond_init(&fetch_condition, NULL);
#endif
	need_fetch_work = false;
}
OOCVolume::~OOCVolume() {
	WakeupFetcher(true);
    if (prefetcher != NULL) {
		int *status;
		prefetcher->join((void**)&status);
    }
}


void OOCVolume::BlockString(const int b[3], char string[1024]) const {
    Volume::BlockString(block_prefix, b, string);
}


bool OOCVolume::Read(const char *fname) {

    gtb::afree<FILE*> f(fopen(fname, "rb"), fclose);

    if (f==0) {
		cerr<<"couldn't open file "<<fname<<endl;
		return false;
    }

    blocksize[0] = 0;
    blocksize[1] = 0;
    blocksize[2] = 0;

    nblocks[0] = 0;
    nblocks[1] = 0;
    nblocks[2] = 0;

    float forigin[3] = {0,0,0};
    translation = Vector3(0,0,0);

    char bp[1024];
    bp[0] = 0;

    dim[0] = -1;
    dim[1] = -1;
    dim[2] = -1;

    float faspect[3] = {0,0,0};
    aspect[0] = 0;
    aspect[1] = 0;
    aspect[2] = 0;

    char line[1024];
    while (fgets(line, 1000, f)) {
		if (3 == sscanf(line, "block size: %d %d %d", &blocksize[0], &blocksize[1], &blocksize[2])) {
		} else if (3 == sscanf(line, "sizes: %d %d %d", &dim[0], &dim[1], &dim[2])) {
		} else if (3 == sscanf(line, "spacings: %g %g %g", &faspect[0], &faspect[1], &faspect[2])) {
			aspect[0]=faspect[0]; aspect[1]=faspect[1]; aspect[2]=faspect[2]; 
		} else if (1 == sscanf(line, "block prefix: %s", bp)) {
		} else if (3 == sscanf(line, "space origin: (%g,%g,%g)", &forigin[0], &forigin[1], &forigin[2])) {
			translation = Vector3(forigin[0], forigin[1], forigin[2]);
		}
    }


    if (dim[0]<0 || dim[1]<0 || dim[2]<0 ||
		aspect[0]==0 || aspect[1]==0 || aspect[2]==0 ||
		bp[0]==0) {
		cerr<<"error reading nooc header"<<endl;
		return false;
    }


    strcpy(block_prefix, fname);
    char *slash = strrchr(block_prefix, '/');

    if (!slash) {
		strcpy(block_prefix, bp);
    } else {
		strcpy(slash+1, bp);
    }


    for (int i=0; i<3; i++) {
		nblocks[i] = dim[i] / blocksize[i];
		if (dim[i] % blocksize[i])
			nblocks[i]++;
    }

    blocks.resize(nblocks[0]);
    for (int i=0; i<nblocks[0]; i++) {
		blocks[i].resize(nblocks[1]);
		for (int j=0; j<nblocks[1]; j++) {
			blocks[i][j].resize(nblocks[2]);
		}
    }

	// int bn=0;
	// for (int i=0; i<nblocks[0]*nblocks[1]*nblocks[2]; i++) {
	// 	cerr<<"block "<<i<<endl;
		
	// 	int b[3];
	// 	while (1) {
	// 		ZCurve::i2xyz(bn++, b[0], b[1], b[2]);
	// 		if (b[0] < nblocks[0] && b[1]<nblocks[1] && b[2]<nblocks[2])
	// 			break;
	// 	}

	// 	cerr<<b[0]<<" "<<b[1]<<" "<<b[2]<<endl;
	// }


    active_index=0;
    active_block[0] = 0;
    active_block[1] = 0;
    active_block[2] = 0;
    active_zindex=0;
    prefetcher = new thlib::Thread(&PrefetchMain_s, (void*)this, 0);
    SetActiveBlock(-1);
    return true;
}

real_type OOCVolume::GetValue(int x, int y, int z) const {

    int block[3] = { x / blocksize[0],
					 y / blocksize[1],
					 z / blocksize[2] };

    int cell[3] = { x % blocksize[0],
					y % blocksize[1],
					z % blocksize[2] };

    // there may be loading/unloading of blocks, so we need to do this atomically
    BlockData &bd = BeginBlockAccess(block);
    real_type ret = bd.rv->GetValue(cell[0], cell[1], cell[2]);
    EndBlockAccess(block);

    return ret;
}


void OOCVolume::Gather(const int cell[3], double nbrs[4][4][4]) const {

	int ncell[3] = { clamp(cell[0],0,dim[0]-1) % blocksize[0],
					 clamp(cell[1],0,dim[1]-1) % blocksize[1],
					 clamp(cell[2],0,dim[2]-1) % blocksize[2] };

	if (ncell[0] > 0 && ncell[0] < blocksize[0]-2 &&
		ncell[1] > 0 && ncell[1] < blocksize[1]-2 &&
		ncell[2] > 0 && ncell[2] < blocksize[2]-2) {
		
		// gather straight out of the block
		int block[3] = { cell[0] / blocksize[0],
						 cell[1] / blocksize[1],
						 cell[2] / blocksize[2] };

		BlockData &bd = BeginBlockAccess(block);
		bd.rv->Gather(ncell, nbrs);
		EndBlockAccess(block);
	} else {

		// otherwise just use the GetValue interface to go to the correct blocks
		for (int z=0; z<4; z++) {
			int zi = clamp(cell[2]+z-1, 0, dim[2]-1);

			for (int y=0; y<4; y++) {
				int yi = clamp(cell[1]+y-1, 0, dim[1]-1);

				for (int x=0; x<4; x++) {
					int xi = clamp(cell[0]+x-1, 0, dim[0]-1);

					nbrs[x][y][z] = (double)GetValue(xi, yi, zi);
				}
			}
		}
    }
}



int OOCVolume::NumBlocks() const {
    return (nblocks[0]*nblocks[1]*nblocks[2]);
}

void OOCVolume::SetActiveBlock(int b) const {

	if (b == active_index) {
		return;
	}

	if (b < 0)
		b = 0;

    int to_increment = b-active_index;
    if (to_increment < 0) {
		to_increment = b;
		active_index = 0;
		active_zindex = 0;
    }

    int64 ai = active_index;
    int64 azi = active_zindex;

    for (int i=0; i<to_increment; i++) {
		ai++;
		while (1) {
			azi++;
			int x,y,z;
			ZCurve::i2xyz(azi, x,y,z);
			if (x<nblocks[0] && y<nblocks[1] && z<nblocks[2])
				break;
		}
    } 


    active_index = ai;
    ZCurve::i2xyz(azi, active_block[0], active_block[1], active_block[2]);
    active_zindex = azi;

	cerr<<"set new active block"<<endl;
	cerr<<"ai: "<<active_index<<endl;
	cerr<<"azi: "<<active_zindex<<endl;
	cerr<<"block: "<<active_block[0]<<" "<<active_block[1]<<" "<<active_block[2]<<endl;

    active_region_loading = true;
	WakeupFetcher();
    while (active_region_loading) { }
}


int OOCVolume::GetActiveBlock() const {
    return active_index;
}

int* OOCVolume::GetActiveBlock3(int i3[3]) const {
	i3[0] = active_block[0];
	i3[1] = active_block[1];
	i3[2] = active_block[2];
	return i3;
}


void OOCVolume::GetBlockBox(int min[3], int max[3]) const {
    for (int i=0; i<3; i++) {
		min[i] = active_block[i]*blocksize[i];
		max[i] = std::min(dim[i], min[i] + blocksize[i]);
    }
}

void OOCVolume::GetBlockBox(int block[3], int min[3], int max[3]) const {
    for (int i=0; i<3; i++) {
		min[i] = block[i]*blocksize[i];
		max[i] = std::min(dim[i], min[i] + blocksize[i]);
    }
}

void OOCVolume::GetBlockMinMax(real_type &min, real_type &max) const {
    BlockData &bd = BeginBlockAccess(active_block);
    bd.rv->GetBlockMinMax(min, max);
    EndBlockAccess(active_block);
}


void OOCVolume::WriteValue(int x, int y, int z, FILE *f) const {

    int block[3] = { x / blocksize[0],
					 y / blocksize[1],
					 z / blocksize[2] };

    int cell[3] = { x % blocksize[0],
					y % blocksize[1],
					z % blocksize[2] };

    // there may be loading/unloading of blocks, so we need to do this atomically
    BlockData &bd = BeginBlockAccess(block);
    bd.rv->WriteValue(cell[0], cell[1], cell[2], f);
    EndBlockAccess(block);
}

#ifdef HAS_ZLIB
void OOCVolume::WriteValue(int x, int y, int z, gzFile f) const {

    int block[3] = { x / blocksize[0],
					 y / blocksize[1],
					 z / blocksize[2] };

    int cell[3] = { x % blocksize[0],
					y % blocksize[1],
					z % blocksize[2] };

    // there may be loading/unloading of blocks, so we need to do this atomically
    BlockData &bd = BeginBlockAccess(block);
    bd.rv->WriteValue(cell[0], cell[1], cell[2], f);
    EndBlockAccess(block);
}
#endif


const char* OOCVolume::TypeString() const {
    int block[3] = {0,0,0};
    BlockData &bd = BeginBlockAccess(block);
    const char *ret = bd.rv->TypeString();
    EndBlockAccess(block);
    return ret;
}


bool OOCVolume::BlockInActiveRegion(const int block[3]) const {

    if (active_index < 0)
		return false;

    extern int prefetch_rings;
    for (int i=0; i<3; i++) {
		int dist = abs(active_block[i]-block[i]);
		if (dist > prefetch_rings)
			return false;
    }

    return true;
}

OOCVolume::BlockData& OOCVolume::BeginBlockAccess(const int block[3]) const {

    BlockData &bd = blocks[block[0]][block[1]][block[2]];

    if (BlockInActiveRegion(block)) {
		if (!bd.rv)
			cerr<<"block in active region not loaded?!?"<<endl;
		return bd;
    } else {

		while (1) {
			cs.enter();
			if (bd.rv) {
				return bd;
			}
			cs.leave();

			for (int i=0; i<3; i++) {
				if (block[i] >= nblocks[i])
					cerr<<"wtfmate"<<endl;
			}

			int64 lb;
			ZCurve::xyz2i(block[0], block[1], block[2], lb);
			blockToLoadASAP = lb;
			cerr<<"loading block ASAP " << lb<<endl;
			//			BREAK;
			WakeupFetcher();
			while (blockToLoadASAP>=0) { }
		}
    }
}


void OOCVolume::EndBlockAccess(const int block[3]) const {
    if (BlockInActiveRegion(block)) {
    } else {
		cs.leave();
    }
}


void OOCVolume::LoadBlock(const int b[3]) {

    char block_string[1024];
    BlockString(b, block_string);
    Volume *v = VolumeRead(block_string);
    if (!v) {
		cerr<<"error reading volume block: "<<block_string<<endl;
		exit(0);
    }

    // make sure the block data matches the ooc header data
    // all we really care about are the dimensions
    int bdim[3];
    BlockDim(b, bdim);
    if (bdim[0] != v->GetDim(0) ||
		bdim[1] != v->GetDim(1) ||
		bdim[2] != v->GetDim(2)) {
		cerr<<"strange block size!"<<endl;
		exit(0);
    }

    if (guidance_computed) {
		vector<Point3> gfpts;
		vector<real_type> gfideals;

		extern char *gf_prefix;
		char fname[1024];
		int ab[3];
		Volume::BlockString(gf_prefix==NULL ? "gf" : gf_prefix,
							b, fname);


		IsoSurfaceGuidanceField::ReadGuidance(fname, gfpts, gfideals);
		v->SetBlockGuidanceField(gfpts, gfideals);
    }


    BlockData &bd = blocks[b[0]][b[1]][b[2]];

    cs.enter();
	if (bd.ordered_traverse_refcount)
		BREAK;
    bd.load_stamp = load_counter++;
    bd.rv = v;
    cs.leave();


	num_block_loads++;
}

void OOCVolume::RemoveLRU(std::set< Array3<int> > &bag) {
    int64 lru_stamp = 0x7fffffffffffffffULL;
    std::set< Array3<int> >::iterator lru_i = bag.end();

    cs.enter();

    for (std::set< Array3<int> >::iterator i=bag.begin(); i!=bag.end(); ++i) {
		const Array3<int> &b = *i;
		const BlockData &bd = blocks[b[0]][b[1]][b[2]];
		if (bd.ordered_traverse_refcount)
			continue;
		int64 i_stamp = bd.load_stamp;
		if (i_stamp < lru_stamp) {
			lru_stamp = i_stamp;
			lru_i = i;
		}
    }

    if (lru_i == bag.end()) {
		cerr<<"no lru to use??"<<endl;
		cs.leave();
		return;
    }


    const Array3<int> &b = *lru_i;
    BlockData &bd_remove = blocks[b[0]][b[1]][b[2]];
    delete bd_remove.rv;
    bd_remove.rv = NULL;
    cs.leave();

    bag.erase(lru_i);
}


bool OOCVolume::FetcherWait() {
#ifdef WIN32
	while (!need_fetch_work) { }
	bool kill = (need_fetch_work == 2);
	need_fetch_work = 0;
	return kill;
#else
	pthread_mutex_lock(&fetch_mutex);

	if (!need_fetch_work) {
		pthread_cond_wait(&fetch_condition, &fetch_mutex);
		if (!need_fetch_work)
			cerr<<"no work for prefetcher?!"<<endl;
	} else {
	}

	bool kill = (need_fetch_work == 2);
	need_fetch_work = 0;
	pthread_mutex_unlock(&fetch_mutex);

	return kill;
#endif
}


void OOCVolume::WakeupFetcher(bool kill) const {
#ifdef WIN32
	need_fetch_work = kill ? 2 : 1;
#else
	pthread_mutex_lock(&fetch_mutex);
	need_fetch_work = kill ? 2 : 1;
	pthread_cond_signal(&fetch_condition);
	pthread_mutex_unlock(&fetch_mutex);
#endif
}


void* OOCVolume::PrefetchMain() {

    extern int max_grab_bag;
    extern int prefetch_rings;
    extern int prefetch_lookahead;

    int prefetched_active = -1;

    std::set< Array3<int> > grab_bag;
    std::set< Array3<int> > prefetched;

    vector< Array3<int> > prefetch_order;
    
    vector< Array3<int> > ring_blocks;

	bool didwork = true;

    while (didwork || !FetcherWait()) {

		didwork = false;

		// for some reason we need to load a block right away
		if (blockToLoadASAP >= 0) {

			Array3<int> block(0,0,0);
			ZCurve::i2xyz(blockToLoadASAP, block[0], block[1], block[2]);
			int ib[3] = { block[0], block[1], block[2] };
			BlockData &bd = blocks[ib[0]][ib[1]][ib[2]];

			if (bd.rv == NULL) {
	    
				bool isgrab = (std::find(prefetch_order.begin(), prefetch_order.end(), block) == prefetch_order.end());

				// if we have to dump something..
				if (isgrab && grab_bag.size()>=max_grab_bag) {
					RemoveLRU(grab_bag);
				}

				if (isgrab) {
					grab_bag.insert(block);
				} else {
					prefetched.insert(block);
					prefetch_order.erase(std::find(prefetch_order.begin(), prefetch_order.end(), block));
				}

 				LoadBlock(ib);
			}

			blockToLoadASAP = -1;
			didwork = true;
		}


		// we've moved to a new active block
		int64 tmpactive = active_zindex;
		if (prefetched_active != tmpactive) {
			prefetched_active = tmpactive;


			// dump everything we have into the grab bag
			grab_bag.insert(prefetched.begin(), prefetched.end());
			prefetched.clear();
			prefetch_order.clear();

			int64 lookahead_index = tmpactive;
			for (int i=0; i<=prefetch_lookahead; i++) {

				int xyz[3];
				ZCurve::i2xyz(lookahead_index, xyz[0],xyz[1],xyz[2]);


				// add this and it's rings to the prefetch list
				for (int ring=0; ring<=prefetch_rings; ring++) {
					GetShellBlocks(xyz, ring, ring_blocks);

					for (unsigned i=0; i<ring_blocks.size(); i++) {
						Array3<int> &b = ring_blocks[i];

						// if we haven't gotten it prefetched yet
						if (prefetched.find(b) == prefetched.end()) {

							// if we've loaded it, but hasn't been marked as prefetch
							if (grab_bag.find(b) != grab_bag.end()) {
								grab_bag.erase(grab_bag.find(b));
								prefetched.insert(b);
							} else {
								// add it to the fetch list if it's not there already
								if (std::find(prefetch_order.begin(), prefetch_order.end(), b) == prefetch_order.end()) {
									prefetch_order.push_back(b);
								}
							}
						}
					}
				}
		

				// go to the next active block we are prefetching
				bool done = false;
				while (1) {
					lookahead_index++;
					int x,y,z;
					ZCurve::i2xyz(lookahead_index, x,y,z);

					// went through all the blocks
					if (x>=nblocks[0] && y>=nblocks[1] && z>=nblocks[2]) {
						done = true;
						break;
					}

					if (x<nblocks[0] && y<nblocks[1] && z<nblocks[2]) {
						break;
					}
				}

				if (done)
					break;
			}

			// empty out the grab bag
			while (grab_bag.size() > max_grab_bag) {
				RemoveLRU(grab_bag);
			}

			didwork=true;
		}


		// if we're waiting for a new active block, check to see if we have
		// everything loaded
		if (active_region_loading) {
			int tmpactive = active_zindex;

			bool all_loaded = true;
			int x,y,z;
			ZCurve::i2xyz(tmpactive, x,y,z);

			for (int xi=std::max(0,x-prefetch_rings); xi<std::min(nblocks[0],x+prefetch_rings+1); xi++) {
				for (int yi=std::max(0,y-prefetch_rings); yi<std::min(nblocks[1],y+prefetch_rings+1); yi++) {
					for (int zi=std::max(0,z-prefetch_rings); zi<std::min(nblocks[2],z+prefetch_rings+1); zi++) {
						Array3<int> b(xi,yi,zi);
						if (prefetched.find(b) == prefetched.end()) {
							all_loaded = false;
						}
					}
				}
			}

			if (all_loaded)
				active_region_loading = false;
		} 


		// finally, load something from the prefetch list
		if (!prefetch_order.empty()) {

			Array3<int> b = prefetch_order[0];
			int ib[3] = { b[0], b[1], b[2] };


			prefetch_order.erase(prefetch_order.begin());
			prefetched.insert(b);

			LoadBlock(ib);

			didwork=true;
		}
    }

    return 0;
}



void OOCVolume::BlockDim(const int b[3], int d[3]) const {
    for (int i=0; i<3; i++) {
		d[i] = blocksize[i];

		// if it's the last block in the dimension
		if (b[i] == nblocks[i]-1) {
			// if the blocks don't divide nicely
			if (blocksize[i] * nblocks[i] != dim[i]) {
				d[i] = dim[i] - (blocksize[i] * (nblocks[i]-1));
			}
		}
    }
}


void DrawBox(const Point3 &min, const Point3 &max, const Point3 &color, bool lines) {


    Box3 box(min, max);
    vector<Point3> pts(5);
    vector<Point3> col(5);

    for (int i=0; i<6; i++) {

		pts.clear();
		col.clear();

		gtb::tPolygon3<real_type> poly = box.face(i);

		for (int j=0; j<4; j++) {
			pts.push_back(poly.point(j));
			col.push_back(color);
		}

		if (lines) {
			pts.push_back(pts.front());
			col.push_back(col.front());
			DbgPLines::add(0, pts, (int)color[0]);
		} else {
			DbgPoly::add(0, pts, col);
		}
    }
}

void OOCVolume::DrawBlock(const Array3<int> &b, const Point3 &color) const {

    Point3 center = Point3((b[0]+0.5)*blocksize[0]*aspect[0],
						   (b[1]+0.5)*blocksize[1]*aspect[1],
						   (b[2]+0.5)*blocksize[2]*aspect[2]);
    center += translation;
    
    DrawBox(center + Vector3(-1,-1,-1),
			center + Vector3( 1, 1, 1),
			color, false);
}

void OOCVolume::ValidateFetchState(std::set<Array3<int> > &prefetched,
								   std::set<Array3<int> > &grab_bag,
								   vector< Array3<int> > &prefetch_list) const {

    for (vector<Array3<int> >::iterator i=prefetch_list.begin(); i!=prefetch_list.end(); ++i) {
		if (prefetched.find(*i) != prefetched.end())
			cerr<<"prefetch list contains already fetched block"<<endl;
		if (grab_bag.find(*i) != grab_bag.end())
			cerr<<"prefetch list contains grab list block"<<endl;
    }

    for (std::set<Array3<int> >::iterator i=prefetched.begin(); i!=prefetched.end(); ++i) {
		if (grab_bag.find(*i) != grab_bag.end())
			cerr<<"grab bag contains prefetched"<<endl;

		if (std::find(prefetch_list.begin(), prefetch_list.end(), *i) != prefetch_list.end())
			cerr<<"prefetch list contains prefetched"<<endl;
    }    

    for (std::set<Array3<int> >::iterator i=grab_bag.begin(); i!= grab_bag.end(); ++i) {
		if (prefetched.find(*i) != prefetched.end())
			cerr<<"prefetched contains grab bag block"<<endl;

		if (std::find(prefetch_list.begin(), prefetch_list.end(), *i) != prefetch_list.end())
			cerr<<"prefetch list contains grab bag block"<<endl;
    }    

}


void OOCVolume::DrawLoaded(std::set<Array3<int> > &prefetched,
						   std::set<Array3<int> > &grab_bag,
						   vector< Array3<int> > &prefetch_list) const {


    dbgClear();


    for (std::set<Array3<int> >::iterator i=prefetched.begin(); i!=prefetched.end(); ++i) {
		int bi[3] = { (*i)[0], (*i)[1], (*i)[2] };
		if (BlockInActiveRegion(bi)) {
			DrawBlock(*i, Point3(1,0,0));
		} else {
			DrawBlock(*i, Point3(0.5,0,0));
		}
    }

    for (std::set<Array3<int> >::iterator i=grab_bag.begin(); i!=grab_bag.end(); ++i) {
		DrawBlock(*i, Point3(0,1,0));
    }

    for (vector<Array3<int> >::iterator i=prefetch_list.begin(); i!=prefetch_list.end(); ++i) {
		DrawBlock(*i, Point3(0,0,1));
    }

    redrawAndWait(' ');
}



void* OOCVolume::OrderedPointTraverseStart(const Point3 &p) {

	OrderedTraverseContext *ctx = new OrderedTraverseContext();
	ctx->p = p;

    real_type dist_to_active = DistToBlock(active_block, p);
    OrderedTraverseInfo oti;
    extern int prefetch_rings;

    if (dist_to_active == 0) {

		// insert the active block and rings
		oti.ring = 0;
		oti.neg_dist_squared = DistToBlock(active_block, p);
		oti.neg_dist_squared = -(oti.neg_dist_squared*oti.neg_dist_squared);
		oti.block = Array3<int>(active_block[0], active_block[1], active_block[2]);
		ctx->pq.push(oti);


		for (int ring=0; ring<=prefetch_rings; ring++) {
			oti.ring = ring+1;
			oti.neg_dist_squared = DistToRingBoundary(p, ring);
			oti.neg_dist_squared = -(oti.neg_dist_squared*oti.neg_dist_squared);
			ctx->pq.push(oti);
		}
    } else {

		// just add all the blocks
		for (int bx=0; bx<nblocks[0]; bx++) {
			for (int by=0; by<nblocks[1]; by++) {
				for (int bz=0; bz<nblocks[2]; bz++) {

					int block[3] = {bx, by, bz};

					oti.ring = 0;
					oti.neg_dist_squared = DistToBlock(block, p);
					oti.neg_dist_squared = -(oti.neg_dist_squared*oti.neg_dist_squared);
					oti.block = Array3<int>(block[0], block[1], block[2]);
					ctx->pq.push(oti);
				}
			}
		}
    }

	return ctx;
}


bool OOCVolume::OrderedPointTraverseNext(void *_ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {

	OrderedTraverseContext *ctx = (OrderedTraverseContext*)_ctx;
	

    while (1) {
	
		if (ctx->pq.empty())
			return false;

		OrderedTraverseInfo top = ctx->pq.top();
		ctx->pq.pop();

		if (top.ring > 0) {

			extern int prefetch_rings;
			if (top.ring > prefetch_rings) {
				cerr<<"not enough rings prefetched for guidance field evaluation!!!"<<endl;
				return false;
			}

			// add all the blocks in this ring
			vector< Array3<int> > shell_blocks;
			GetShellBlocks(active_block, top.ring, shell_blocks);

			for (unsigned i=0; i<shell_blocks.size(); i++) {
				OrderedTraverseInfo oti;
				oti.ring = 0;
				oti.neg_dist_squared = DistToBlock(shell_blocks[i].x, ctx->p);
				oti.neg_dist_squared = -(oti.neg_dist_squared*oti.neg_dist_squared);
				oti.block = shell_blocks[i];
				ctx->pq.push(oti);
			}
		} else {

			std::map<Array3<int>,void*>::iterator si = ctx->subctx.find(top.block);
			void *subctx = (si == ctx->subctx.end()) ? NULL : (*si).second;

			// see if we have to start a new traverse for the block
			BlockData &bd = BeginBlockAccess(top.block.x);

			bool started = (subctx != NULL);
			if (!started) {
				cs.enter();
				bd.ordered_traverse_refcount++;
				cs.leave();
				subctx = ctx->subctx[top.block] = bd.rv->OrderedPointTraverseStart(ctx->p);
			}

			// set up the next point from this block
			Point3 b_p;
			real_type b_distsq, b_ideal;

			bool blockfinished = !bd.rv->OrderedPointTraverseNext(subctx, b_distsq, b_p, b_ideal);
			if (!blockfinished) {
				OrderedTraverseInfo oti;
				oti.ring = 0;
				oti.neg_dist_squared = -b_distsq;
				oti.block = top.block;
				oti.p = b_p;
				oti.ideal = b_ideal;
				ctx->pq.push(oti);
			}
			EndBlockAccess(top.block.x);

			if (started) {
				// we can actually return what we got
				squared_dist = -top.neg_dist_squared;
				point = top.p;
				ideal = top.ideal;
				return true;
			}
		}
    }
}


void OOCVolume::OrderedPointTraverseEnd(void *_ctx) {

	OrderedTraverseContext *ctx = (OrderedTraverseContext*)_ctx;

    for (std::map<Array3<int>,void*>::iterator i=ctx->subctx.begin(); i != ctx->subctx.end(); ++i) {
		BlockData &bd = BeginBlockAccess((*i).first.x);
		bd.rv->OrderedPointTraverseEnd((*i).second);
		cs.enter();
		bd.ordered_traverse_refcount--;
		cs.leave();
		EndBlockAccess((*i).first.x);
    }

	delete ctx;
}


void OOCVolume::SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal) {
    BlockData &bd = BeginBlockAccess(active_block);

    bd.rv->SetBlockGuidanceField(pts, ideal);

    EndBlockAccess(active_block);

	// only write it out if we're not reusing it
	if (!reuse_guidance) {
		extern char *gf_prefix;
		char fname[1024];
		Volume::BlockString(gf_prefix==NULL ? "gf" : gf_prefix,
							active_block, fname);
		IsoSurfaceGuidanceField::WriteGuidance(fname, pts, ideal);
	}
}

void OOCVolume::SetGuidanceFieldDone() {

	// make sure we have the gf points for all loaded block loaded!
	for (int bx=0; bx<blocks.size(); bx++) {
		for (int by=0; by<blocks[bx].size(); by++) {
			for (int bz=0; bz<blocks[bx][by].size(); bz++) {

				BlockData &bd = blocks[bx][by][bz];
				if (bd.rv != NULL) {

					vector<Point3> gfpts;
					vector<real_type> gfideals;

					extern char *gf_prefix;
					char fname[1024];
					int b[3] = {bx,by,bz};
					Volume::BlockString(gf_prefix==NULL ? "gf" : gf_prefix,
										b, fname);

					IsoSurfaceGuidanceField::ReadGuidance(fname, gfpts, gfideals);
					bd.rv->SetBlockGuidanceField(gfpts, gfideals);
				}
			}
		}
	}

    guidance_computed = true;
}


void OOCVolume::SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds) {

	extern char *gf_prefix;
	char seed_prefix[1024];
    char fname[1024];
	sprintf(seed_prefix, "%s.seeds", (gf_prefix==NULL ? "gf" : gf_prefix));
    Volume::BlockString(seed_prefix, active_block, fname);

	if (!reuse_seeds) {
		WriteSeeds(fname, bnd_pts, bnd_norms, loops, cc_seeds);
	}
}

void OOCVolume::GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const {

	extern char *gf_prefix;
	char seed_prefix[1024];
    char fname[1024];
	sprintf(seed_prefix, "%s.seeds", (gf_prefix==NULL ? "gf" : gf_prefix));
    Volume::BlockString(seed_prefix, active_block, fname);

	ReadSeeds(fname, bnd_pts, bnd_norms, loops, cc_seeds);
}


void OOCVolume::GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
    SetActiveBlock(0);

    int bmins[3], bmaxs[3];
    GetBlockBox(bmins, bmaxs);

    box = Box3(Point3(1e34,1e34,1e34), Point3(-1e34,-1e34,-1e34));
    for (int i=0; i<3; i++) {
		if (GetAspect(i) > 0) {
			box.set_min(i, (bmins[i]-1)*GetAspect(i) + translation[i]);
			box.set_max(i, (bmaxs[i]+1)*GetAspect(i) + translation[i]);
		} else {
			box.set_max(i, (bmins[i]-1)*GetAspect(i) + translation[i]);
			box.set_min(i, (bmaxs[i]+1)*GetAspect(i) + translation[i]);
		}
    }


	GetBlockSeeds(pts, norms, loops, seeds);
}


bool OOCVolume::GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
    int ab = GetActiveBlock();
    if (ab == NumBlocks()-1)
		return false;
    SetActiveBlock(ab+1);

    int bmins[3], bmaxs[3];
    GetBlockBox(bmins, bmaxs);
    box = Box3(Point3(1e34,1e34,1e34), Point3(-1e34,-1e34,-1e34));
    for (int i=0; i<3; i++) {
		if (GetAspect(i) > 0) {
			box.set_min(i, (bmins[i]-1)*GetAspect(i) + translation[i]);
			box.set_max(i, (bmaxs[i]-1)*GetAspect(i) + translation[i]);
		} else {
			box.set_max(i, (bmins[i]-1)*GetAspect(i) + translation[i]);
			box.set_min(i, (bmaxs[i]-1)*GetAspect(i) + translation[i]);
		}
    }


	GetBlockSeeds(pts, norms, loops, seeds);

    return true;
}

int OOCVolume::GetBlockForPoint(const Point3 &p) const {

	int cell[3];

    Point3 pt = p - translation;
    for (int i=0; i<3; i++) {
		cell[i] = round_to_minus_inf(pt[i]/aspect[i]);
    }

	//	double xl[3];	// local coordinate
	//	LocalInfo(p, cell, xl);

	int block[3];
	for (int i=0; i<3; i++) {
		block[i] = cell[i] / blocksize[i];
	}

	return std::max(0, ((block[2] * nblocks[1]) + block[1]) * nblocks[0] + block[0]);
}



void OOCVolume::GetShellBlocks(int block[3], int shell, vector< Array3<int> > &shell_blocks) const {

    shell_blocks.clear();

    int zstart = std::max(0,          block[2]-shell);
    int zend   = std::min(nblocks[2], block[2]+shell+1);

    int ystart = std::max(0, block[1]-shell);
    int yend   = std::min(nblocks[1], block[1]+shell+1);

    int xstart = std::max(0, block[0]-shell);
    int xend   = std::min(nblocks[0], block[0]+shell+1);

    for (int z=zstart; z<zend; z++) {
		bool zshell = (z==block[2]-shell || z==block[2]+shell);
		for (int y=ystart; y<yend; y++) {
			bool yshell = (y==block[1]-shell || y==block[1]+shell);
			for (int x=xstart; x<xend; x++) {
				bool xshell = (x==block[0]-shell || x==block[0]+shell);

				if (zshell || yshell || xshell){
					shell_blocks.push_back(Array3<int>(x,y,z));
				}
			}
		}
    }
}


real_type OOCVolume::DistToRingBoundary(const Point3 &p, int ring) const {

    real_type dist = 1e34;

    for (int i=0; i<3; i++) {

		int bmin = active_block[i] - ring;
		int bmax = active_block[i] + ring + 1;

		real_type fmin = ((bmin>0)          ? (bmin*blocksize[i]) : -1e34) * aspect[i] + translation[i];
		real_type fmax = ((bmax<nblocks[i]) ? (bmax*blocksize[i]) :  1e34) * aspect[i] + translation[i];
		if (fmin > fmax)
			std::swap(fmin,fmax);

		if (p[i] < fmin) return 0;
		if (p[i] > fmax) return 0;

		real_type d = std::min(gtb::abs(p[i]-fmin), gtb::abs(p[i]-fmax));
		dist = std::min(dist, d);
    }

    return dist;
}


real_type OOCVolume::DistToBlock(int block[3], const Point3 &p) const {
    int bmin[3], bmax[3];
    GetBlockBox(block, bmin, bmax);

    real_type fmin[3], fmax[3];
    for (int i=0; i<3; i++) {
		if (aspect[i] > 0) {
			fmin[i] = bmin[i]*aspect[i] + translation[i];
			fmax[i] = bmax[i]*aspect[i] + translation[i];
		} else {
			fmin[i] = bmax[i]*aspect[i] + translation[i];
			fmax[i] = bmin[i]*aspect[i] + translation[i];
		}
    }

    Box3 box(Point3(fmin[0],fmin[1],fmin[2]), Point3(fmax[0],fmax[1],fmax[2]));
    //    cerr<<"contains"<<box.contains(p)<<endl;
    //    cerr<<"distance"<<box.distance(p)<<endl;
    return box.distance(p);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Volume* VolumeRead(const char *fname) {

    bool isregular=false;
    bool isooc = false;

    NHDR nhdr;
    nhdr.sizes[0] = nhdr.sizes[1] = DIM;
    nhdr.sizes[2] = 0;
    nhdr.spacings[0] = nhdr.spacings[1] = 1;
    nhdr.spacings[2] = 2;
    nhdr.endian = false;
    nhdr.translation = Vector3(0,0,0);
    strcpy(nhdr.datafile, fname);

    // read / create a nhdr for the given file
    bool endswith(const char *s, const char *end);
    if (endswith(fname, ".short.vol")) {
		nhdr.type = NHDR::Short;
		isregular=true;
    } else if (endswith(fname, ".byte.vol")) {
		nhdr.type = NHDR::UChar;
		isregular=true;
    } else if (endswith(fname, ".float.vol")) {
		nhdr.type = NHDR::Float;
		isregular=true;
    } else if (endswith(fname, ".nhdr")) {
		if (!ReadNHDR(fname, nhdr))
			return NULL;
		isregular=true;
    } else if (endswith(fname, ".nooc")) {
		isooc = true;
    } else {
		cerr<<"unknown volume format: "<<fname<<endl;
		return NULL;
    }


    // read the actual data
    if (isregular) {
		switch (nhdr.type) {
		case NHDR::Short: {
			RegularVolume<short> *v = new RegularVolume<short>();
			if (!v->ReadNRRD(nhdr)) {
				delete v;
				return NULL;
			}
			return v;
		}
		case NHDR::Float: {
			RegularVolume<float> *v = new RegularVolume<float>();
			if (!v->ReadNRRD(nhdr)) {
				delete v;
				return NULL;
			}
			return v;
		}
		case NHDR::UChar: {
			RegularVolume<unsigned char> *v = new RegularVolume<unsigned char>();
			if (!v->ReadNRRD(nhdr)) {
				delete v;
				return NULL;
			}
			return v;
		}
		}
    } else if (isooc) {

		OOCVolume *v = new OOCVolume();
		if (!v->Read(fname)) {
			delete v;
			return NULL;
		}
		return v;
    }

    return NULL;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////



IsoSurfaceProjector::IsoSurfaceProjector(const Volume &vol, real_type iso) : volume(vol), isovalue(iso) {

}


IsoSurfaceProjector::~IsoSurfaceProjector() {
}


int IsoSurfaceProjector::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const {
    return ProjectPoint(fp, tp, tn, NULL);
}


int IsoSurfaceProjector::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn, const Vector3 *plane) const {

    // newton step onto the surface
    tp = fp;
    int last_cell[3]={-1,-1,-1};
    double nbrs[4][4][4];

    extern int max_iter;
    extern real_type value_tol;
    extern real_type step_tol;

    real_type max_dist_sq = Vector3(volume.GetAspect(0),
									volume.GetAspect(1),
									volume.GetAspect(2)).squared_length();

    int iter=0;
    while (1) {

		if (Point3::squared_distance(fp,tp) > max_dist_sq*4) {
			return PROJECT_FAILURE;
		}

		int this_cell[3];
		double xl[3];	// local coordinate
		volume.LocalInfo(tp, this_cell, xl);


		if (this_cell[0]!=last_cell[0] || this_cell[1]!=last_cell[1] || this_cell[2]!=last_cell[2]) {

			volume.Gather(this_cell, nbrs);
			last_cell[0]=this_cell[0];
			last_cell[1]=this_cell[1];
			last_cell[2]=this_cell[2];
		}


		double gradient[3];
		volume.spline.Gradient(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), nbrs, xl, gradient);
		double f = volume.spline.Eval(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), nbrs, xl);

		tn[0]=(real_type)gradient[0];
		tn[1]=(real_type)gradient[1];
		tn[2]=(real_type)gradient[2];


		if (plane) {
			double dot = gradient[0]*(*plane)[0] + gradient[1]*(*plane)[1] + gradient[2]*(*plane)[2];
			gradient[0] -= dot*(*plane)[0];
			gradient[1] -= dot*(*plane)[1];
			gradient[2] -= dot*(*plane)[2];
		}


		double gradmag = (gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2]);
		if (gradmag == 0) {
			cerr<<" g";
			return PROJECT_FAILURE;
		}
		
		double step = -(f-isovalue) / (gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2]);

		tp[0] += (real_type)(step*gradient[0]);
		tp[1] += (real_type)(step*gradient[1]);
		tp[2] += (real_type)(step*gradient[2]);

		// check for stopping
		bool smallstep = step<step_tol;
		bool smallval  = fabs(f-isovalue)<value_tol;
	
		if (smallstep && smallval) {
			break;
		}

		iter++;
		if (iter>max_iter) {
			cerr<<" ";
			cerr<< (smallstep ? "s" : "S");
			cerr<< (smallval  ? "v" : "V");
			return PROJECT_FAILURE;
		}

    }

    tn.normalize();
    return PROJECT_SUCCESS;
}



int IsoSurfaceProjector::ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn, const Point3 &center, const Vector3 &axis) const {

	real_type radius = Point3::distance(fp, center);

    Vector3 udir = fp - center;
    udir.normalize();
    Vector3 vdir = axis.cross(udir);
    vdir.normalize();

    udir *= radius;
    vdir *= radius;

    RingEval re(*this, center, udir, vdir);


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

			if (theta_max>M_PI_2) {
				//		cerr<<"couldn't bracket positive"<<endl;
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
      
      
			if (theta_min<-M_PI_2) {
				//		cerr<<"couldn't bracket negative"<<endl;
				return PROJECT_FAILURE;
			}
      
		} while (f_theta_min>0);
    
    }

  
    extern real_type step_tol;
    real_type zero = NR::zbrent(re, theta_min, theta_max, step_tol*radius);
    tp = re.point(zero);


    int this_cell[3];
    double nbrs[4][4][4];
    double xl[3];	// local coordinate
    if (!volume.LocalInfo(tp, this_cell, xl)) return PROJECT_FAILURE; // BOUNDARY
    volume.Gather(this_cell, nbrs);

    double gradient[3];
    volume.spline.Gradient(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), nbrs, xl, gradient);
    tn[0]=(real_type)gradient[0];
    tn[1]=(real_type)gradient[1];
    tn[2]=(real_type)gradient[2];
    tn.normalize();

	if (tn[0]!=tn[0] ||
		tn[1]!=tn[1] ||
		tn[2]!=tn[2])
		return PROJECT_FAILURE;


    return PROJECT_SUCCESS;
}



int IsoSurfaceProjector::ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const {

#ifdef ISO_CIRCULAR_PROJECTION

    Vector3 axis = base2.position - base1.position;
    axis.normalize();
  
    Point3 center = base1.position + (axis.dot(fp - base1.position)) * axis;

	int res = ProjectPoint(fp, tp, tn, center, axis);

#else

    int res = ProjectPoint(fp, tp, tn);

#endif


	if (res == PROJECT_SUCCESS && curvature!=NULL) {

		real_type k1, k2;
		if (volume.CurvatureAtPoint(tp, k1, k2)) {
			*curvature = k2/3;
			if (fabs(k1)>fabs(k2))
				*curvature = k1/3;
		}
	}

	return res;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool gfbound_twomethod = false; // seems faster this way
void IsoSurfaceGuidanceField::SampleCell(int id, int cell[3], 
										 double nbrs[4][4][4], double coefs[4][4][4],
										 const Interval xspan_s[3], IsoSurfaceProjector &projector,
										 vector<Point3> &pts, vector<real_type> &ideals,
										 const int _third_sign[3][3][3],
										 int depth) const {


	// check for quick termination
	Interval fspan = Interval::Intersect(volume.spline.Eval2(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), coefs, xspan_s),
										  volume.spline.Eval(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), nbrs, xspan_s));
	if (fspan.max<projector.GetIsoValue() || fspan.min>projector.GetIsoValue())
		return; // no intersection at all!


	// third partial tensor
	int third_sign[3][3][3];
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			for (int k=0; k<3; k++) {
				third_sign[i][j][k] = _third_sign[i][j][k];
			}
		} 
	}

	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {
			for (int k=j; k<3; k++) {

				if (_third_sign[i][j][k] == 0) {

					int partial[3] = {0,0,0};
					partial[i]++;
					partial[j]++;
					partial[k]++;

					Interval third_partial = volume.spline.EvalPartial2(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), coefs, xspan_s, partial);

					if (third_partial.max < 0) {
						third_sign[i][j][k] = -1;
						third_sign[i][k][j] = -1;
						third_sign[j][i][k] = -1;
						third_sign[j][k][i] = -1;
						third_sign[k][i][j] = -1;
						third_sign[k][j][i] = -1;
					}

					if (third_partial.min > 0) {
						third_sign[i][j][k] = 1;
						third_sign[i][k][j] = 1;
						third_sign[j][i][k] = 1;
						third_sign[j][k][i] = 1;
						third_sign[k][i][j] = 1;
						third_sign[k][j][i] = 1;
					}
				}
			}
		}
	}


	// use 3rd partial tensor to put tighter bounds on the hessian
	gtb::tmat3<Interval> hessian;
	for (int i=0; i<3; i++) {
		for (int j=i; j<3; j++) {

			int partial[3] = { 0, 0, 0 };
			partial[i]++;
			partial[j]++;

			Interval nspan = volume.spline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
												 nbrs, coefs, partial, xspan_s,
												 third_sign[i][j][0], third_sign[i][j][1], third_sign[i][j][2], gfbound_twomethod);
			hessian[i][j] = nspan;
			hessian[j][i] = nspan;
		}
	}
 

	// use hessian to put tighter bounds on the gradient
	Interval gradient[3];
	for (int i=0; i<3; i++) {
		int partial[3] = { 0, 0, 0 };
		partial[i] = 1;
		gradient[i] = volume.spline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
										  nbrs, coefs, partial, xspan_s,
										  interval_sign(hessian[i][0]), interval_sign(hessian[i][1]), interval_sign(hessian[i][2]), gfbound_twomethod);
	}


	// use gradient to put tighter bounds on f
	{
		int partial[3] = { 0, 0, 0 };
		fspan = volume.spline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
									nbrs, coefs, partial, xspan_s,
									interval_sign(gradient[0]), interval_sign(gradient[1]), interval_sign(gradient[2]), gfbound_twomethod);
	}

	if (fspan.max<projector.GetIsoValue() || fspan.min>projector.GetIsoValue())
		return; // no intersection at all!



	// do ||A^k||^(1/k)
	// with k=1 gives some benefit - going too far will overflow to inf's
	int k = 1;
	gtb::tmat3<Interval> h_k = hessian;
	for (int p=0; p<k; p++) {
		gtb::tmat3<Interval> h_k1 = h_k;
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				h_k[i][j] = h_k1[i][0]*h_k1[0][j] + h_k1[i][1]*h_k1[1][j] + h_k1[i][2]*h_k1[2][j];
			}
		}
	}

 
	// try 2 different norms and pick the better one
	Interval frobnorm = sqrt(square(h_k[0][0]) + 2*square(h_k[0][1]) + 2*square(h_k[0][2]) +
							 square(h_k[1][1]) + 2*square(h_k[1][2]) + square(h_k[2][2]));

	Interval rownorm = Interval::Union(abs(h_k[0][0]) + abs(h_k[0][1]) + abs(h_k[0][2]),
									   Interval::Union(abs(h_k[1][0]) + abs(h_k[1][1]) + abs(h_k[1][2]),
													   abs(h_k[2][0]) + abs(h_k[2][1]) + abs(h_k[2][2])));
									   

	for (int p=0; p<k; p++) {
		frobnorm = sqrt(frobnorm);
		rownorm = sqrt(rownorm);
	}

	Interval norm = Interval::Intersect(frobnorm, rownorm);



	Interval gradmag = sqrt(square(gradient[0]) + square(gradient[1]) + square(gradient[2]));

	real_type diag = sqrt(square((xspan_s[0][1]-xspan_s[0][0])) + 
						  square((xspan_s[1][1]-xspan_s[1][0])) +
						  square((xspan_s[2][1]-xspan_s[2][0])));

	// paper has 2*sqrt(3) that is not needed
	//real_type kmax_bound = ((2*sqrt(3.f) * norm / gradmag).max);
	real_type kmax_bound = ((norm / gradmag).max);
	real_type ideal_bound = MaxCurvatureToIdeal(kmax_bound/3);


	// subdivide if needed
	int samples_needed = (gradmag.min<=0) ? 1000000 :  std::max(1,((int)(diag/ideal_bound / 2)));
	if (samples_needed > curvature_sub) {
		real_type center[3] = { (xspan_s[0].min+xspan_s[0].max)/2,
								(xspan_s[1].min+xspan_s[1].max)/2,
								(xspan_s[2].min+xspan_s[2].max)/2 };

		for (int x=0; x<2; x++) {
			real_type xmin = (x==0) ? xspan_s[0].min : center[0];
			real_type xmax = (x==0) ? center[0] : xspan_s[0].max;

			for (int y=0; y<2; y++) {
				real_type ymin = (y==0) ? xspan_s[1].min : center[1];
				real_type ymax = (y==0) ? center[1] : xspan_s[1].max;
				
				for (int z=0; z<2; z++) {
					real_type zmin = (z==0) ? xspan_s[2].min : center[2];
					real_type zmax = (z==0) ? center[2] : xspan_s[2].max;

					Interval subi[3] = { Interval(xmin,xmax),
										 Interval(ymin,ymax),
										 Interval(zmin,zmax) };

					SampleCell(id, cell, nbrs, coefs, subi, projector, pts, ideals, third_sign, depth+1);
				}
			}
		}

	} else {


 		real_type spacing = 1.0 / samples_needed;

		Point3 min(xspan_s[0].min + cell[0] * volume.GetAspect(0),
				   xspan_s[1].min + cell[1] * volume.GetAspect(1),
				   xspan_s[2].min + cell[2] * volume.GetAspect(2));
		Point3 max(xspan_s[0].max + cell[0] * volume.GetAspect(0),
				   xspan_s[1].max + cell[1] * volume.GetAspect(1),
				   xspan_s[2].max + cell[2] * volume.GetAspect(2));
		min += volume.GetTranslation();
		max += volume.GetTranslation();
				   
 
 		for (int x=0; x<samples_needed; x++) {
 			for (int y=0; y<samples_needed; y++) {
 				for (int z=0; z<samples_needed; z++) {

 					Point3 sl((x+0.5)*spacing,
 							  (y+0.5)*spacing,
 							  (z+0.5)*spacing);

					sl[0] = sl[0]*(max[0]-min[0]) + min[0];
					sl[1] = sl[1]*(max[1]-min[1]) + min[1];
					sl[2] = sl[2]*(max[2]-min[2]) + min[2];


					// If we *REALLY* wanted a 'sufficient' gf, we would add (sl,kmax_bound), but that adds a TON of points
					// and is extra conservative (slooooow and small triangles).
					// pts.push_back(sl);
					// ideals.push_back(MaxCurvatureToIdeal(kmax_bound/3));
					
 					Point3 tp;
 					Vector3 tn;
 					real_type k1,k2;
					
 					if (projector.ProjectPoint(sl, tp, tn)==PROJECT_SUCCESS) {
						if (volume.CurvatureAtPoint(tp, k1, k2)) {

							real_type kmax = k2;
							if (fabs(k1)>fabs(k2))
								kmax = k1;

							if (tp[0]>=min[0] && tp[0]<=max[0] &&
								tp[1]>=min[1] && tp[1]<=max[1] &&
								tp[2]>=min[2] && tp[2]<=max[2]) {

								pts.push_back(tp);
								ideals.push_back(MaxCurvatureToIdeal(kmax/3));
							} else {
								// We don't care about points that project outside of the interval we're looking at.
							}

						}
 					}
 				}
 			}
 		}
	}
}



void IsoSurfaceGuidanceField::SampleCell(int id, int cell[3], vector<a444_wrap> &nbrs, vector<a444_wrap> &coefs, int considerbits,
										 const Interval xspan_s[3], IsoSurfaceProjector &projector, MultiMaterialVolume &volume,
										 vector<Point3> &pts, vector<real_type> &ideals, 
										 int depth) const {

	vector<Interval> fspans(nbrs.size());

	for (int n=0; n<nbrs.size(); n++) {
		if (!(considerbits&(1<<n)))
			continue;

		fspans[n] = Interval::Intersect(volume.cspline.Eval2(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), coefs[n].n, xspan_s),
										volume.cspline.Eval(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), nbrs[n].n, xspan_s));



		// third
		Interval third[3][3][3];
		for (int i=0; i<3; i++) {
			for (int j=i; j<3; j++) {
				for (int k=j; k<3; k++) {

					int partial[3] = { 0, 0, 0 };
					partial[i]++;
					partial[j]++;
					partial[k]++;

					Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
														  nbrs[n].n, coefs[n].n, partial, xspan_s,
														  0,0,0,
														  true);

					third[i][j][k] = nspan;
					third[i][k][j] = nspan;
					third[j][i][k] = nspan;
					third[j][k][i] = nspan;
					third[k][i][j] = nspan;
					third[k][j][i] = nspan;
				}
			}
		}



		// hessian
		gtb::tmat3<Interval> hessian;
		for (int i=0; i<3; i++) {
			for (int j=i; j<3; j++) {

				int partial[3] = { 0, 0, 0 };
				partial[i]++;
				partial[j]++;

				Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
													  nbrs[n].n, coefs[n].n, partial, xspan_s,
													  interval_sign(third[i][j][0]),
													  interval_sign(third[i][j][1]),
													  interval_sign(third[i][j][2]),
													  true);
				hessian[i][j] = nspan;
				hessian[j][i] = nspan;
			}
		}



		// gradient
		Interval grad[3];
		for (int i=0; i<3; i++) {
			int partial[3] = { 0, 0, 0 };
			partial[i]++;

			Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
												  nbrs[n].n, coefs[n].n, partial, xspan_s,
												  interval_sign(hessian[i][0]),
												  interval_sign(hessian[i][1]),
												  interval_sign(hessian[i][2]),
												  true);
			grad[i] = nspan;
		}


		// value
		int partial[3] = { 0, 0, 0 };
		fspans[n] = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2),
										 nbrs[n].n, coefs[n].n, partial, xspan_s,
										 interval_sign(grad[0]),
										 interval_sign(grad[1]),
										 interval_sign(grad[2]),
										 true);
	}


	real_type max_min = -INFINITY;
	real_type min_max = INFINITY;
	for (int i=0; i<nbrs.size(); i++) {
		if (!(considerbits&(1<<i)))
			continue;

		max_min = std::max(max_min, fspans[i].min);
		min_max = std::min(min_max, fspans[i].max);
	}


	int myconsiderbits = 0;
	for (int i=0; i<nbrs.size(); i++) {
		if (!(considerbits&(1<<i)))
			continue;

		if (fspans[i].max >= max_min)
			myconsiderbits |= 1<<i;
	}


	int nconsiders = 0;
	vector<int> considers;
	for (int i=0; i<nbrs.size(); i++) {
		if (myconsiderbits&(1<<i)) {
			nconsiders++;
			considers.push_back(i);
		}
	}



	// look more closely at the pairings
	if (nconsiders > 2) {
		for (int i=0; i<considers.size()-1; i++) {
			for (int j=i+1; j<considers.size(); j++) {

				int ci = considers[i];
				int cj = considers[j];

				real_type tnbrs[4][4][4];
				real_type tcoefs[4][4][4];

				for (int k=0; k<4*4*4; k++) {
					(&tnbrs[0][0][0])[k] = (&nbrs[ci].n[0][0][0])[k] - (&nbrs[cj].n[0][0][0])[k];
					(&tcoefs[0][0][0])[k] = (&coefs[ci].n[0][0][0])[k] - (&coefs[cj].n[0][0][0])[k];
				}



				// third
				Interval third[3][3][3];
				for (int i=0; i<3; i++) {
					for (int j=i; j<3; j++) {
						for (int k=j; k<3; k++) {

							int partial[3] = { 0, 0, 0 };
							partial[i]++;
							partial[j]++;
							partial[k]++;

							Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
																  tnbrs, tcoefs, partial, xspan_s,
																  0,0,0,
																  true);

							third[i][j][k] = nspan;
							third[i][k][j] = nspan;
							third[j][i][k] = nspan;
							third[j][k][i] = nspan;
							third[k][i][j] = nspan;
							third[k][j][i] = nspan;
						}
					}
				}


				// hessian
				gtb::tmat3<Interval> hessian;
				for (int i=0; i<3; i++) {
					for (int j=i; j<3; j++) {

						int partial[3] = { 0, 0, 0 };
						partial[i]++;
						partial[j]++;

						Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
															  tnbrs, tcoefs, partial, xspan_s,
															  interval_sign(third[i][j][0]),
															  interval_sign(third[i][j][1]),
															  interval_sign(third[i][j][2]),
															  true);
						hessian[i][j] = nspan;
						hessian[j][i] = nspan;
					}
				}



				// gradient
				Interval grad[3];
				for (int i=0; i<3; i++) {
					int partial[3] = { 0, 0, 0 };
					partial[i]++;

					Interval nspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2), 
														  tnbrs, tcoefs, partial, xspan_s,
														  interval_sign(hessian[i][0]),
														  interval_sign(hessian[i][1]),
														  interval_sign(hessian[i][2]),
														  true);
					grad[i] = nspan;
				}


				// value
				int partial[3] = { 0, 0, 0 };
				Interval tspan = volume.cspline.Bound(volume.GetAspect(0), volume.GetAspect(1), volume.GetAspect(2),
													  tnbrs, tcoefs, partial, xspan_s,
													  interval_sign(grad[0]),
													  interval_sign(grad[1]),
													  interval_sign(grad[2]),
													  true);


				if (tspan.min > 0) {
					// i - j > 0  ->  j is never in this sub-cell
					//					cerr<<j<<" never in cell"<<endl;

					myconsiderbits &= ~(1<<cj);
				}
			
				if (tspan.max < 0) {
					// i - j < 0  ->  i is never in this sub-cell
					//					cerr<<i<<" never in cell"<<endl;

					myconsiderbits &= ~(1<<ci);
				}
			
			}
		}
	}



	nconsiders = 0;
	considers.clear();
	for (int i=0; i<nbrs.size(); i++) {
		if (myconsiderbits&(1<<i)) {
			nconsiders++;
			considers.push_back(i);
		}
	}	 


	if (nconsiders <= 1) // less than 2 materials in this sub-cell
		return;

	// if (nconsiders <= 2) // at most 2 materials in this sub-cell
	// 	return;

	// if we go too deep into an intersection, bail out
	if (depth >= 2 && nconsiders != 2)
		return;


	if (depth >= 2 && nconsiders == 2) {
		Point3 center = Point3(volume.GetAspect(0) * cell[0] + 0.5*(xspan_s[0][0]+xspan_s[0][1]),
							   volume.GetAspect(1) * cell[1] + 0.5*(xspan_s[1][0]+xspan_s[1][1]),
							   volume.GetAspect(2) * cell[2] + 0.5*(xspan_s[2][0]+xspan_s[2][1])) + volume.GetTranslation();

		Point3 fp = center;

		volume.cur_mat1 = considers[0];
		volume.cur_mat2 = considers[1];

		// try projecting the point onto all the surfaces we think run through here
		Point3 tp;
		Vector3 tn;
		if (projector.ProjectPoint(fp, tp, tn) != PROJECT_SUCCESS) {
			return;
		}

		fp = tp;


		real_type k1,k2;
		if (volume.CurvatureAtPoint(fp, k1, k2)) {
			
			real_type kmax = k2;
			if (fabs(k1)>fabs(k2))
				kmax = k1;

			pts.push_back(fp);
			ideals.push_back(MaxCurvatureToIdeal(kmax/3));
		}

		return;
	}
	


	// got here we need to subdivide
	real_type center[3] = { (xspan_s[0].min+xspan_s[0].max)/2,
							(xspan_s[1].min+xspan_s[1].max)/2,
							(xspan_s[2].min+xspan_s[2].max)/2 };

	for (int x=0; x<2; x++) {
		real_type xmin = (x==0) ? xspan_s[0].min : center[0];
		real_type xmax = (x==0) ? center[0] : xspan_s[0].max;

		for (int y=0; y<2; y++) {
			real_type ymin = (y==0) ? xspan_s[1].min : center[1];
			real_type ymax = (y==0) ? center[1] : xspan_s[1].max;
				
			for (int z=0; z<2; z++) {
				real_type zmin = (z==0) ? xspan_s[2].min : center[2];
				real_type zmax = (z==0) ? center[2] : xspan_s[2].max;

				Interval subi[3] = { Interval(xmin,xmax),
									 Interval(ymin,ymax),
									 Interval(zmin,zmax) };

				SampleCell(id, cell, nbrs, coefs, myconsiderbits, subi, projector, volume, pts, ideals, depth+1);
			}
		}
	}
}





bool IsoSurfaceGuidanceField::AtMostSingleIntersection(double nbrs[4][4][4], const Point3 &x1, const Point3 &x2, real_type isovalue) const {

	bool xeq = (x1[0]==x2[0]);
	bool yeq = (x1[1]==x2[1]);
	bool zeq = (x1[2]==x2[2]);

	if (!xeq && yeq && zeq) {

		// along x axis
		Interval x(x1[0], x2[0]);
		real_type y = x1[1];
		real_type z = x1[2];

		// check if the function value crosses the isovalue along the edge
		int partial[3] = {0,0,0};
		Interval f = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);

		if (f.max < isovalue || f.min > isovalue)
			return true; // don't span

		// we may have an intersection - check if the gradient along the edge is always increasing or decreasing
		partial[0]++;
		Interval g = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);

		if (g.max <= 0 || g.min >= 0)
			return true; // monotonic along edge

	} else if (xeq && !yeq && zeq) {

		// along y axis
		real_type x = x1[0];
		Interval y(x1[1], x2[1]);
		real_type z = x1[2];

		// check if the function value crosses the isovalue along the edge
		int partial[3] = {0,0,0};
		Interval f = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);

		if (f.max < isovalue || f.min > isovalue)
			return true; // don't span

		// we may have an intersection - check if the gradient along the edge is always increasing or decreasing
		partial[1]++;
		Interval g = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);

		if (g.max <= 0 || g.min >= 0)
			return true; // monotonic along edge

	} else if (xeq && yeq && !zeq) {

		// along z axis
		real_type x = x1[0];
		real_type y = x1[1];
		Interval z(x1[2], x2[2]);

		// check if the function value crosses the isovalue along the edge
		int partial[3] = {0,0,0};
		Interval f = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);
		if (f.max < isovalue || f.min > isovalue)
			return true; // don't span

		// we may have an intersection - check if the gradient along the edge is always increasing or decreasing
		partial[2]++;
		Interval g = volume.spline.EvalPartial(1, 1, 1, nbrs, x, y, z, partial);
		if (g.max <= 0 || g.min >= 0)
			return true; // monotonic along edge

	} else {
		cerr<<"AtMostSingleIntersection() not axis aligned????"<<endl;
		exit(1);
	}


	return false;
}


int maxdepth = 0;
int IsoSurfaceGuidanceField::BoundaryRecursionDepth(int axis, double nbrs[4][4][4], double coefs[4][4][4], real_type isovalue, 
													const Interval xspan[3], int depth) const {


	if (depth>maxdepth) {
		cerr<<"new max depth: "<<depth<<endl;
 		maxdepth = depth;
	}

	// if (depth > 10)
	// 	return 10;


	if (xspan[axis][0] != xspan[axis][1]) {
		cerr<<"wtf xspan not collapsed?"<<endl;
	}

	int axis1 = (axis+1)%3;
	int axis2 = (axis+2)%3;


	// gradient
	int grad_signs[3];

	for (int i=0; i<3; i++) {

		if (i==axis) {
			grad_signs[i] = 1;

		} else {

			int partial[3]= {0,0,0};
			partial[i]++;

			grad_signs[i] = interval_sign(volume.spline.Bound(1,1,1, nbrs, coefs, partial, xspan, 0,0,0, true));
		}
	}



	// check for f going into the face
	{
		int partial[3] = { 0, 0, 0 };
		Interval fspan = volume.spline.Bound(1,1,1, nbrs, coefs, partial, xspan, grad_signs[0], grad_signs[1], grad_signs[2], true);

		if (fspan.max<isovalue || fspan.min>isovalue)
			return 0; // no intersection at all - don't even need this level
	}



	bool must_recurse = false;

	// if all of the gradients change sign, we need to recurse
	if (grad_signs[axis1]==0 && grad_signs[axis2]==0) {
		must_recurse = true;


	} else {

		// otherwise one of the partials is 0 across the face - so the boundary must just be a line (no circles)
		// make sure each edge only has a single crossing at most
		

		if (!must_recurse) {
			Point3 x1, x2;
			x1[axis] = x2[axis] = xspan[axis][0];

			x1[axis1] = xspan[axis1][0];
			x2[axis1] = xspan[axis1][0];

			x1[axis2] = xspan[axis2][0];
			x2[axis2] = xspan[axis2][1];

			if (!AtMostSingleIntersection(nbrs, x1, x2, isovalue))
				must_recurse = true;
		}

		if (!must_recurse) {
			Point3 x1, x2;
			x1[axis] = x2[axis] = xspan[axis][0];

			x1[axis1] = xspan[axis1][1];
			x2[axis1] = xspan[axis1][1];

			x1[axis2] = xspan[axis2][0];
			x2[axis2] = xspan[axis2][1];

			if (!AtMostSingleIntersection(nbrs, x1, x2, isovalue))
				must_recurse = true;
		}

		if (!must_recurse) {
			Point3 x1, x2;
			x1[axis] = x2[axis] = xspan[axis][0];

			x1[axis1] = xspan[axis1][0];
			x2[axis1] = xspan[axis1][1];

			x1[axis2] = xspan[axis2][0];
			x2[axis2] = xspan[axis2][0];

			if (!AtMostSingleIntersection(nbrs, x1, x2, isovalue))
				must_recurse = true;
		}

		if (!must_recurse) {
			Point3 x1, x2;
			x1[axis] = x2[axis] = xspan[axis][0];

			x1[axis1] = xspan[axis1][0];
			x2[axis1] = xspan[axis1][1];

			x1[axis2] = xspan[axis2][1];
			x2[axis2] = xspan[axis2][1];

			if (!AtMostSingleIntersection(nbrs, x1, x2, isovalue))
				must_recurse = true;
		}

	}


	// if we don't need to recurse - then this level is fine
	if (!must_recurse)
		return depth+1;

	// otherwise - subdivide
	real_type center[3] = { (xspan[0].min+xspan[0].max)/2,
							(xspan[1].min+xspan[1].max)/2,
							(xspan[2].min+xspan[2].max)/2 };


	int maxdepth = depth;
	for (int a1=0; a1<2; a1++) {
		real_type a1min = (a1==0) ? xspan[axis1].min : center[axis1];
		real_type a1max = (a1==0) ? center[axis1] : xspan[axis1].max;

		for (int a2=0; a2<2; a2++) {
			real_type a2min = (a2==0) ? xspan[axis2].min : center[axis2];
			real_type a2max = (a2==0) ? center[axis2] : xspan[axis2].max;

			Interval subi[3];
			subi[axis] = xspan[axis];
			subi[axis1] = Interval(a1min, a1max);
			subi[axis2] = Interval(a2min, a2max);

			maxdepth = std::max(maxdepth, BoundaryRecursionDepth(axis, nbrs, coefs, isovalue, subi, depth+1));
		}
	}

	return maxdepth;
}



class LineEvalCell {

    public:
    LineEvalCell(const TrivariateSpline<double> &_spline, const double (&_nbrs)[4][4][4], const Point3 &_p0, const Point3 &_p1, real_type iso) :
		spline(_spline), nbrs(_nbrs), p0(_p0), p1(_p1), iso_value(iso) { }

    real_type operator()(real_type t) const {

		Point3 x = p0 + t*(p1-p0);
	    
		double e = spline.Eval(1,1,1, nbrs, &x[0]);

		return (real_type)(e - iso_value);
    }


	const TrivariateSpline<double> &spline;
    const double (&nbrs)[4][4][4];
    const real_type iso_value;
    const Point3 p0;
    const Point3 p1;
};


void IsoSurfaceGuidanceField::GetBlockFaceBoundaries(real_type isovalue, int axis, bool min, vector< vector<Point3> > &bsegments, vector<bool> &bloops) const { 

	int bmin[3];
	int bmax[3];
	volume.GetBlockBox(bmin, bmax);

	int axis1 = (axis+1)%3;
	int axis2 = (axis+2)%3;


	std::map<int,int> edges;
	vector<int> edges2;

	std::map<int,Point3> sharedpts;
	vector<Point3> unsharedpts;

 	vector<int> veryboundpts;


	int maxaxis1 = std::min(bmax[axis1], volume.GetDim(axis1)-1);
	int maxaxis2 = std::min(bmax[axis2], volume.GetDim(axis2)-1);

	for (int i=bmin[axis1]; i<maxaxis1; i++) {
		for (int j=bmin[axis2]; j<maxaxis2; j++) {
			
			int cell[3] = { 0,0,0 };
			cell[axis] = min ? bmin[axis] : std::min(bmax[axis], volume.GetDim(axis)-1);
			cell[axis1] = i;
			cell[axis2] = j;


			// get the depth we need to go to to be accurate
			double nbrs[4][4][4];
			volume.Gather(cell, nbrs);

			double coefs[4][4][4];
			volume.spline.Precompute(nbrs, coefs);


			Interval xspan[3] = { Interval(0,1), Interval(0,1), Interval(0,1) };
			xspan[axis][1] = xspan[axis][0];


			int depth = BoundaryRecursionDepth(axis, nbrs, coefs, isovalue, xspan, 0);

			if (depth == 0)
				continue;

			int divisions = pow(2,depth-1);


			std::map<int, int> indexmap; // for unshared indices

			// create the vertices we need
			int nsharedi1=0;
			int nsharedi2=0;
			int nsharedj1=0;
			int nsharedj2=0;
			for (int i2=0; i2<=divisions; i2++) {
				for (int j2=0; j2<=divisions; j2++) {

					// add along axis1
					if (i2 != divisions) {

						real_type l1[3];
						l1[axis]  = 0;
						l1[axis1] = (i2/(real_type)(divisions));
						l1[axis2] = (j2/(real_type)(divisions));
						real_type fl1 = volume.spline.Eval(1,1,1, nbrs, l1) - isovalue;

						real_type l2[3];
						l2[axis]  = 0;
						l2[axis1] = ((i2+1)/(real_type)(divisions));
						l2[axis2] = (j2/(real_type)(divisions));
						real_type fl2 = volume.spline.Eval(1,1,1, nbrs, l2) - isovalue;


						if (fl1*fl2 < 0) {
							
							LineEvalCell le(volume.spline, nbrs, l1, l2, isovalue);

							extern real_type step_tol;
							extern real_type value_tol;
							real_type step = step_tol;

							real_type s;
							for (int t=0; t<5; t++) {
								s = NR::zbrent(le, (real_type)0, (real_type)1, step);
								step *= 0.05;
								if (fabs(le(s)) < value_tol)
									break;
							}


							if (j2==0 || j2==divisions) {
								// shared point
								
								if (j2==0 || j==maxaxis2-1) {

									int index;
									if (j2 == 0) {
										if (nsharedi1==3) {
											cerr<<"more than 3 shared verts on edge???"<<endl;
											exit(1);
										}

										index = (j*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedi1;
										nsharedi1++;
									} else {
										if (nsharedi2==3) {
											cerr<<"more than 3 shared verts on edge???"<<endl;
											exit(1);
										}

										index = ((j+1)*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedi2;
										nsharedi2++;
									}


									Point3 p = Point3::inbetween(Point3(l2), Point3(l1), s) + Vector3(cell[0], cell[1], cell[2]);
									p[0] *= volume.GetAspect(0);
									p[1] *= volume.GetAspect(1);
									p[2] *= volume.GetAspect(2);
									p += volume.GetTranslation();

									sharedpts[index] = p;
									indexmap[(j2*(divisions+1) + i2)*2 + 0] = index; // positive

									if (j==bmin[axis2] && j2==0  || j==maxaxis2-1 && j2==divisions)
										veryboundpts.push_back(index);
								}

							} else {
								// unshared point
								Point3 p = Point3::inbetween(Point3(l2), Point3(l1), s) + Vector3(cell[0], cell[1], cell[2]);
								p[0] *= volume.GetAspect(0);
								p[1] *= volume.GetAspect(1);
								p[2] *= volume.GetAspect(2);
								p += volume.GetTranslation();

								indexmap[(j2*(divisions+1) + i2)*2 + 0] = -unsharedpts.size()-1; // negative
								unsharedpts.push_back(p);
							}
						}

					}


					// add along axis2
					if (j2 != divisions) {

						real_type l1[3];
						l1[axis]  = 0;
						l1[axis1] = (i2/(real_type)(divisions));
						l1[axis2] = (j2/(real_type)(divisions));
						real_type fl1 = volume.spline.Eval(1,1,1, nbrs, l1) - isovalue;

						real_type l2[3];
						l2[axis]  = 0;
						l2[axis1] = (i2/(real_type)(divisions));
						l2[axis2] = ((j2+1)/(real_type)(divisions));
						real_type fl2 = volume.spline.Eval(1,1,1, nbrs, l2) - isovalue;


						if (fl1*fl2 < 0) {
							
							LineEvalCell le(volume.spline, nbrs, l1, l2, isovalue);

							extern real_type step_tol;
							extern real_type value_tol;
							real_type step = step_tol;

							real_type s;
							for (int t=0; t<5; t++) {
								s = NR::zbrent(le, (real_type)0, (real_type)1, step);
								step *= 0.05;
								if (fabs(le(s)) < value_tol)
									break;
							}


							if (i2==0 || i2==divisions) {
								// shared point
								
								if (i2==0 || i==maxaxis1-1) {

									int index;
									if (i2 == 0) {
										if (nsharedj1==3) {
											cerr<<"more than 3 shared verts on edge???"<<endl;
											exit(1);
										}

										index = (j*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedj1 + 3;
										nsharedj1++;
									} else {
										if (nsharedj2==3) {
											cerr<<"more than 3 shared verts on edge???"<<endl;
											exit(1);
										}

										index = (j*(maxaxis1-bmin[axis1]+1) + (i+1))*6 + nsharedj2 + 3;
										nsharedj2++;
									}


									Point3 p = Point3::inbetween(Point3(l2), Point3(l1), s) + Vector3(cell[0], cell[1], cell[2]);
									p[0] *= volume.GetAspect(0);
									p[1] *= volume.GetAspect(1);
									p[2] *= volume.GetAspect(2);
									p += volume.GetTranslation();

									sharedpts[index] = p;
									indexmap[(j2*(divisions+1) + i2)*2 + 1] = index; // positive

									if (i==bmin[axis1] && i2==0  || i==maxaxis1-1 && i2==divisions)
										veryboundpts.push_back(index);
								}

							} else {
								// unshared point
								Point3 p = Point3::inbetween(Point3(l2), Point3(l1), s) + Vector3(cell[0], cell[1], cell[2]);
								p[0] *= volume.GetAspect(0);
								p[1] *= volume.GetAspect(1);
								p[2] *= volume.GetAspect(2);
								p += volume.GetTranslation();

								indexmap[(j2*(divisions+1) + i2)*2 + 1] = -unsharedpts.size()-1; // negative
								unsharedpts.push_back(p);
							}
						}

					}

				}
			}


			// create the connections
			nsharedi1=0;
			nsharedi2=0;
			nsharedj1=0;
			nsharedj2=0;
			for (int i2=0; i2<divisions; i2++) {
				for (int j2=0; j2<divisions; j2++) {

					real_type imin = (i2/(real_type)(divisions));
					real_type imax = ((i2+1)/(real_type)(divisions));

					real_type jmin = (j2/(real_type)(divisions));
					real_type jmax = ((j2+1)/(real_type)(divisions));


					int vindexes[4];

					if (j2==0) {
						vindexes[0] = (j*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedi1;
					} else {
						vindexes[0] = indexmap[(j2*(divisions+1) + i2)*2 + 0];
					}

					if (i2==(divisions-1)) {
						vindexes[1] = (j*(maxaxis1-bmin[axis1]+1) + (i+1))*6 + nsharedj2 + 3;
					} else {
						vindexes[1] = indexmap[(j2*(divisions+1) + (i2+1))*2 + 1];
					}

					if (j2==(divisions-1)) {
						vindexes[2] = ((j+1)*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedi2;
					} else {
						vindexes[2] = indexmap[((j2+1)*(divisions+1) + i2)*2 + 0];
					}

					if (i2==0) {
						vindexes[3] = (j*(maxaxis1-bmin[axis1]+1) + i)*6 + nsharedj1 + 3;
					} else {
						vindexes[3] = indexmap[(j2*(divisions+1) + i2)*2 + 1];
					}

					


 
					Point3 quad[4];
					for (int k=0; k<4; k++) {
						quad[k][axis] = 0;
					}

					quad[0][axis1] = imin;
					quad[0][axis2] = jmin;

					quad[1][axis1] = imax;
					quad[1][axis2] = jmin;

					quad[2][axis1] = imax;
					quad[2][axis2] = jmax;

					quad[3][axis1] = imin;
					quad[3][axis2] = jmax;


					real_type fquad[4];
					for (int k=0; k<4; k++) {
						real_type localx[3] = { quad[k][0],
												quad[k][1],
												quad[k][2] };
						fquad[k] = volume.spline.Eval(1,1,1, nbrs, localx);
					}
					

					for (int k=0; k<4; k++) {
						if (fquad[k]<isovalue && fquad[(k+1)%4]>isovalue) {
							// entering here - look for where we're leaving

							if (k==0 && j2==0) {
								nsharedi1++;
							}
							if (k==1 && i2==(divisions-1)) {
								nsharedj2++;
							}
							if (k==2 && j2==(divisions-1)) {
								nsharedi2++;
							}
							if (k==3 && i2==0) {
								nsharedj1++;
							}


							for (int _k2=1; _k2<4; _k2++) {
								int k2 = (k+_k2)%4;

								if (fquad[k2]>isovalue && fquad[(k2+1)%4]<isovalue) {
									// leaving here!


									if (k2==0 && j2==0) {
										nsharedi1++;
									}
									if (k2==1 && i2==(divisions-1)) {
										nsharedj2++;
									}
									if (k2==2 && j2==(divisions-1)) {
										nsharedi2++;
									}
									if (k2==3 && i2==0) {
										nsharedj1++;
									}


									edges2.push_back(vindexes[k]);
									edges2.push_back(vindexes[k2]);

									edges[vindexes[k]] = vindexes[k2];

									//									cerr<<vindexes[k]<<" "<<vindexes[k2]<<endl;

									break;
								}

							}

						}
					}

				}
			}

		}
	}


	// stitch the edges together
	vector< vector<int> > ibounds;
	vector<bool> loops;

	// start with all the points that exit the boundary
	for (int i=0; i<veryboundpts.size(); i++) {
		if (edges.find(veryboundpts[i]) == edges.end())
			continue;

		vector<int> e;

		e.push_back(veryboundpts[i]);

		while (1) {

			if (edges.find(e.back()) == edges.end())
				break; // no edges going off the last vert

			int next = edges[e.back()];

			edges.erase(e.back());

			e.push_back(next);
		}

		ibounds.push_back(e);
		loops.push_back(false);
	}

	while (edges.size()) {

		vector<int> e;

		e.push_back(edges.begin()->first);

		while (1) {

			if (edges.find(e.back()) == edges.end())
				break; // no edges going off the last vert

			int next = edges[e.back()];

			edges.erase(e.back());

			e.push_back(next);
		}

		if (e.front() != e.back()) {
			cerr<<"loop not a loop??"<<endl;
		}

		e.pop_back();

		ibounds.push_back(e);
		loops.push_back(true);
	}


	for (int b=0; b<ibounds.size(); b++) {

		vector<Point3> seg(ibounds[b].size());

		for (int i=0; i<ibounds[b].size(); i++) {
			if (ibounds[b][i] < 0) {
				seg[i] = unsharedpts[-ibounds[b][i]-1]; // unshared
			} else {
				seg[i] = sharedpts[ibounds[b][i]]; // shared
			}
		}


		if (!min)
			reverse(seg);

		bsegments.push_back(seg);
		bloops.push_back(loops[b]);
	}
}



void IsoSurfaceGuidanceField::BuildGuidanceParallel(int nt, int id, CSObject &cs, real_type isovalue, IsoSurfaceProjector &projector, vector<Point3> &pts, vector<real_type> &ideals) const { 

    vector<Point3> mypoints;
    vector<real_type> myideals;

	mypoints.reserve(10000);
	myideals.reserve(10000);

	int bmin[3];
	int bmax[3];
	volume.GetBlockBox(bmin, bmax);

	for (int z=bmin[2]+id; z<bmax[2]; z+=nt) {
		for (int y=bmin[1]; y<bmax[1]; y++) {
			for (int x=bmin[0]; x<bmax[0]; x++) {

				if (x == volume.GetDim(0)-1)    continue;
				if (y == volume.GetDim(1)-1)    continue;
				if (z == volume.GetDim(2)-1)    continue;

				int cell[3] = {x,y,z};
				double nbrs[4][4][4];
				volume.Gather(cell, nbrs);
				Interval cspan[3] = { Interval(0,1) * volume.GetAspect(0),
									  Interval(0,1) * volume.GetAspect(1),
									  Interval(0,1) * volume.GetAspect(2) };

				double coefs[4][4][4];
				volume.spline.Precompute(nbrs, coefs);


				int third_sign[3][3][3];
				{
					for (int i=0; i<3; i++) {
						for (int j=0; j<3; j++) {
							for (int k=0; k<3; k++) {
								third_sign[i][j][k] = 0;
							}
						}
					}
				}

				SampleCell(id, cell, nbrs, coefs, cspan, projector, mypoints, myideals, third_sign);
			}
		}
	}


    cs.enter();
	pts.insert(pts.end(), mypoints.begin(), mypoints.end());
	ideals.insert(ideals.end(), myideals.begin(), myideals.end());
    cs.leave();
}



void IsoSurfaceGuidanceField::BuildMultiGuidanceParallel(int nt, int id, CSObject &cs, MultiMaterialVolume &volume, IsoSurfaceProjector &projector, vector<Point3> &pts, vector<real_type> &ideals) const { 

    vector<Point3> mypoints;
    vector<real_type> myideals;

    mypoints.reserve(10000);
    myideals.reserve(10000);

	int bmin[3];
	int bmax[3];
	volume.GetBlockBox(bmin, bmax);

	for (int z=bmin[2]+id; z<bmax[2]; z+=nt) {
		for (int y=bmin[1]; y<bmax[1]; y++) {
			for (int x=bmin[0]; x<bmax[0]; x++) {

				if (x == volume.GetDim(0)-1)    continue;
				if (y == volume.GetDim(1)-1)    continue;
				if (z == volume.GetDim(2)-1)    continue;

				int cell[3] = {x,y,z};

				Interval cspan[3] = { Interval(0,1) * volume.GetAspect(0),
				 					  Interval(0,1) * volume.GetAspect(1),
									  Interval(0,1) * volume.GetAspect(2) };


				int considerbits = 0;
				vector<a444_wrap> nbrs;
				vector<a444_wrap> coefs;
				for (int i=0; i<volume.matProbs.size(); i++) {
					nbrs.push_back(a444_wrap());
					coefs.push_back(a444_wrap());

					volume.matProbs[i]->Gather(cell, nbrs.back().n);
					volume.spline.Precompute(nbrs.back().n, coefs.back().n);

					considerbits |= 1<<i;
				}


				SampleCell(id, cell, nbrs, coefs, considerbits, cspan, projector, volume, mypoints, myideals);
			}
		}
	}

    cs.enter();
    for (unsigned i=0; i<mypoints.size(); i++)
		pts.push_back(mypoints[i]);
    for (unsigned i=0; i<myideals.size(); i++)
		ideals.push_back(myideals[i]);
    cs.leave();
}



void SplitBoundaryLists(const vector<int> &bound, const vector<int> &types, 
						vector< vector<int> > &bsegments, vector<bool> &loops, vector<int> &seg_types) {

    for (int t=0; t<6; t++) {
		int type = 1<<t;

		for (unsigned i=0; i<bound.size(); i++) {
			if (types[bound[i]] == 0)
				cerr<<"non-type!"<<endl;
			unsigned p = (i+bound.size()-1) % bound.size();

			if ((types[bound[i]] & type) && !(types[bound[p]] & type)) {

				// this is the start of a segment on type's boundary
				vector<int> seg;
				for (unsigned j=0; j<bound.size(); j++) {
					unsigned n = (i+j)%bound.size();
					if (!(types[bound[n]] & type))
						break;
					seg.push_back(bound[n]);
				}

		
				if (seg.size() == bound.size())
					cerr<<"wtfmate"<<endl;

				loops.push_back(false);

				seg_types.push_back(type);
				bsegments.push_back(seg);
			}
		}
    }

    // check for a single loop
    if (loops.size() == 0) {
		loops.push_back(true);
		seg_types.push_back(types[bound[0]]);
		bsegments.push_back(bound);
    }
}


void GetConnectedComponents(const TriangleMesh &mesh, vector< vector<int> > &components) {

    vector<bool> flooded(mesh.verts.size());
    for (unsigned i=0;i<mesh.verts.size(); i++) {
		flooded[i]=false;
    }

    while (1) {
	    
		int startv=-1;
		for (unsigned i=0; i<flooded.size(); i++) {
			if (!flooded[i]) {
				if (mesh.verts[i].someface < 0) {
					flooded[i] = true;
				} else {
					startv=i; break;
				}
			}
		}
		if (startv<0) break;

		if (mesh.verts[startv].someface < 0)
			cerr<<"someface < 0"<<endl;


		vector<int> component_points;
		component_points.push_back(startv);

	    
		vector<int> toflood;
		toflood.push_back(startv);
		flooded[startv] = true;
	    
		while (!toflood.empty()) {
			int tf = toflood.back();
			toflood.pop_back();
			for (TriangleMesh::VertexVertexIteratorI vi(mesh, tf); !vi.done(); ++vi) {
				int tf2 = *vi;
				if (flooded[tf2]) continue;
				flooded[tf2]=true;
				component_points.push_back(tf2);
				toflood.push_back(tf2);
			}
		}

		components.push_back(component_points);
    }
}

void ProjectListToIntersection(const IsoSurfaceProjector &volproj,
							   vector<Point3> &pts, vector<Vector3> &norms,
							   const Vector3 &plane, real_type cell_diag, bool loop) {

	
    extern real_type value_tol;
    vector<Point3> r_pts;
    vector<Vector3> r_norms;

    if (loop) {
		Point3 p = pts[0];
		Vector3 n = norms[0];
		pts.push_back(p);
		norms.push_back(n);
    }

	// {
	// 	dbgClear();
	// 	for (unsigned i=0; i<pts.size(); i++) {
	// 		if (loop) {
	// 			DbgPoints::add(pts[i], 0,1,0);
	// 		} else {
	// 			DbgPoints::add(pts[i], 1,0,0);
	// 		}

	// 	}
	// 	if (loop)
	// 		DbgPLines::add(pts,1);
	// 	else
	// 		DbgPLines::add(pts,0);

	// 	redrawAndWait(' ');
	// }

    for (unsigned i=0; i<pts.size()-1; i++) {

		r_pts.push_back(pts[i]);
		r_norms.push_back(norms[i]);


		Point3 start = pts[i];
		Point3 end = pts[i+1];

		int num = (int)(40 * Point3::distance(start, end) / cell_diag);

		for (int j=1; j<num; j++) {
			real_type x = (real_type)j / num;

			Point3 fp = start + x*(end-start);
			Point3 tp;
			Vector3 tn;

			if (volproj.ProjectPoint(fp, tp, tn, start, plane) == PROJECT_SUCCESS) {
				r_pts.push_back(tp);
				r_norms.push_back(tn);
			}
		}
    }

    r_pts.push_back(pts.back());
    r_norms.push_back(norms.back());


	// for (unsigned i=0; i<r_pts.size(); i++) {
	// 	DbgPoints::add(r_pts[i], 0,0,1);
	// }
	// DbgPLines::add(r_pts, 2);
	// redrawAndWait(' ');



    // put the points back, if they didn't project to some place crazy
    pts.clear();
    norms.clear();


    pts.push_back(r_pts[0]);
    norms.push_back(r_norms[0]);
    for (unsigned i=1; i<r_pts.size()-1; i++) {
		real_type idist = std::max(Point3::distance(r_pts[i-1], r_pts[i]),
								   Point3::distance(r_pts[i+1], r_pts[i]));
		real_type ndist = Point3::distance(r_pts[i-1], r_pts[i+1]);

		if (idist < ndist) {
			pts.push_back(r_pts[i]);
			norms.push_back(r_norms[i]);
		}
	    

    }
    pts.push_back(r_pts.back());
    norms.push_back(r_norms.back());


    if (loop) {
		pts.pop_back();
		norms.pop_back();
    }
}



IsoSurfaceGuidanceField::IsoSurfaceGuidanceField(IsoSurfaceProjector &projector, MultiMaterialVolume &vol, vector< vector<Point3> > &loops, const vector< vector<int> > &loop_mats,
												 real_type rho, real_type min_step, real_type max_step, real_type reduction)
    : GuidanceField(rho, min_step, max_step, reduction), volume(vol) {

	vector<Point3> gf_pts;
	vector<real_type> gf_ideals;

	vector<Point3> tpts;
	vector<real_type> tideals;

	// densely sample the loops
	const int density = 10;
	for (int l=0; l<loops.size(); l++) {

		cerr<<"loop "<<l<<endl;
		for (int i=0; i<loop_mats[l].size(); i++)
			cerr<<loop_mats[l][i]<<" ";
		cerr<<endl;

		if (loops[l].size() == 1) {
			cerr<<"cc seed"<<endl;
			continue;
		}


		vector<Point3> nloop;
		nloop.reserve(loops[l].size() * density);

		for (int i=0; i<loops[l].size(); i++) {
			Point3 start = loops[l][i];
			Point3 end = loops[l][(i+1)%loops[l].size()];

			for (int j=0; j<density; j++) {
				Point3 fp = Point3::inbetween(start, end, 1 - (real_type)j / density);

				for (int n=0; n<4; n++) {
					for (int i=0; i<loop_mats[l].size()-1; i++) {
						for (int j=i+1; j<loop_mats[l].size(); j++) {

							int ci = loop_mats[l][i];
							int cj = loop_mats[l][j];

							vol.cur_mat1 = ci;
							vol.cur_mat2 = cj;

							Point3 tp;
							Vector3 tn;
							if (projector.ProjectPoint(fp, tp, tn) != PROJECT_SUCCESS) {
								goto skip_loop_point;
							}

							fp = tp;
						}
					}		
				}

				nloop.push_back(fp);

				skip_loop_point:
				continue;
			}

		}

		// put the dense loop back so we can resample it for initial boundaries
		loops[l] = nloop;

	}


	if (!reuse_guidance) {
		cerr<<"projecting gf points"<<endl;
		ParallelExecutor(idealNumThreads, makeClassFunctor(this,&IsoSurfaceGuidanceField::BuildMultiGuidanceParallel), vol, projector, gf_pts, gf_ideals);
		int num_orig = gf_pts.size();


		// compute curvatures of loops
		for (int l=0; l<loops.size(); l++) {
			if (loops[l].size() == 1)
				continue;

			vector<Point3> &nloop = loops[l];
			vector<real_type> ideals(nloop.size());

			for (int j=0; j<nloop.size(); j++) {

				int min, max;
				min = j-5;
				max = j+5;

				vector<Point3> epts;
				vector<real_type> weights;

				for (int k=min; k<=max; k++) {
					epts.push_back(nloop[(k+nloop.size())%nloop.size()]);
					weights.push_back(1);
				}

				real_type k = edge_curvature(epts, weights, j-min);
				gf_pts.push_back(nloop[j]);
				gf_ideals.push_back(MaxCurvatureToIdeal(k));
			}
		}


		// trim 1
		vector<Point3> tpts1;
		vector<real_type> tideals1;

		vector<int> marked1(gf_ideals.size());
		Trimmer< GetPointVector<Point3> > trimmer1(*this, GetPointVector<Point3>(gf_pts), gf_ideals);
		trimmer1.Trim(marked1);

		for (unsigned i=0; i<marked1.size(); i++) {
			if (!marked1[i]) {
				tpts1.push_back(gf_pts[i]);
				tideals1.push_back(gf_ideals[i]);
			}
		}

		// trim 2
		vector<int> marked(tideals1.size());
		Trimmer< GetPointVector<Point3> > trimmer(*this, GetPointVector<Point3>(tpts1), tideals1);
		trimmer.Trim(marked);

		for (unsigned i=0; i<marked.size(); i++) {
			if (!marked[i]) {
				tpts.push_back(tpts1[i]);
				tideals.push_back(tideals1[i]);
			}
		}

		WriteGuidance("tmp.gf", tpts, tideals);
	} else {
		ReadGuidance("tmp.gf", tpts, tideals);
	}


	volume.SetBlockGuidanceField(tpts, tideals);


	// resample all the loops
	for (int l=0; l<loops.size(); l++) {
		if (loops[l].size() == 1)
			continue;

		vector< vector<Vector3> > in, on;
		vector<Point3> op;
		ResampleCurve(loops[l], in, op, on, false);
		loops[l] = op;
	}
	
}



IsoSurfaceGuidanceField::IsoSurfaceGuidanceField(IsoSurfaceProjector &projector, Volume &vol, real_type rho, real_type min_step, real_type max_step, real_type reduction, int _block)
    : GuidanceField(rho, min_step, max_step, reduction), volume(vol) {

    int num_orig = 0;
    int num_trimmed=0;

    real_type cell_diag = Vector3(volume.GetAspect(0), 
								  volume.GetAspect(1), 
								  volume.GetAspect(2)).length();



    real_type isovalue = projector.GetIsoValue();

    for (int block=(_block<0 ? 0 : _block);
		 block< (_block<0 ? volume.NumBlocks() : _block+1);
		 block++) {

		volume.SetActiveBlock(block);

		cerr<<"computing guidance field for block "<<block<<" of "<<volume.NumBlocks()<<endl;
		int ab[3];
		volume.GetActiveBlock3(ab);




		if (!reuse_seeds) {

			int block_min[3];
			int block_max[3];
			volume.GetBlockBox(block_min, block_max);

			int block_type = 0;
			if (block_min[0] == 0)                   block_type |= MC_VERT_TYPE_XMIN;
			if (block_max[0] == volume.GetDim(0))    block_type |= MC_VERT_TYPE_XMAX;
			if (block_min[1] == 0)                   block_type |= MC_VERT_TYPE_YMIN;
			if (block_max[1] == volume.GetDim(1))    block_type |= MC_VERT_TYPE_YMAX;
			if (block_min[2] == 0)                   block_type |= MC_VERT_TYPE_ZMIN;
			if (block_max[2] == volume.GetDim(2))    block_type |= MC_VERT_TYPE_ZMAX;


			cerr<<"finding boundaries"<<endl;
			vector< vector<Point3> > bsegments;
			vector<bool> bloops;
			vector<int> types;

			if (block_type & MC_VERT_TYPE_XMIN) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 0, true, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_XMIN);
			}

			if (block_type & MC_VERT_TYPE_YMIN) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 1, true, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_YMIN);
			}

			if (block_type & MC_VERT_TYPE_ZMIN) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 2, true, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_ZMIN);
			}

			if (block_type & MC_VERT_TYPE_XMAX) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 0, false, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_XMAX);
			}

			if (block_type & MC_VERT_TYPE_YMAX) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 1, false, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_YMAX);
			}

			if (block_type & MC_VERT_TYPE_ZMAX) {
				GetBlockFaceBoundaries(projector.GetIsoValue(), 2, false, bsegments, bloops);
				while (types.size() < bloops.size())
					types.push_back(MC_VERT_TYPE_ZMAX);
			}


			vector< vector<Point3> > used_points;
			vector< vector<Vector3> > used_norms;
			vector<bool> used_loops;


			for (unsigned i=0; i<bsegments.size(); i++) {

				vector<Point3> seg;
				vector<Vector3> norms;

				for (unsigned j=0; j<bsegments[i].size(); j++) {
					seg.push_back(bsegments[i][j]);

					Vector3 normal;
					volume.NormalAtPoint(seg.back(), normal);
					norms.push_back(normal);
				}

				Vector3 plane(0,0,0);

				if (types[i] & MC_VERT_TYPE_XMIN)
					plane[0] = -1;
				if (types[i] & MC_VERT_TYPE_XMAX)
					plane[0] = 1;
				if (types[i] & MC_VERT_TYPE_YMIN)
					plane[1] = -1;
				if (types[i] & MC_VERT_TYPE_YMAX)
					plane[1] = 1;
				if (types[i] & MC_VERT_TYPE_ZMIN)
					plane[2] = -1;
				if (types[i] & MC_VERT_TYPE_ZMAX)
					plane[2] = 1;

				if (plane.length() != 1)
					cerr<<"not axis aligned plane!"<<endl;

				extern bool allow_outside;
				allow_outside=true;
				int size;
				do {
					size = seg.size();
					ProjectListToIntersection(projector, seg, norms, plane, cell_diag, bloops[i]);
				} while (size < seg.size());
				allow_outside=false;


				if (seg.size() > 1) {
					used_points.push_back(seg);
					used_norms.push_back(norms);
					used_loops.push_back(bloops[i]);
				}
			}


			// run marching cubes to find connected components
			cerr<<"marching cubes"<<endl;
			TriangleMesh mc_mesh;
			vector<int> vertex_types;
			{
				OutputControllerITS oc_its;
				ControllerWrapper cw(NULL, NULL, &oc_its);	// to catch the mc output
				MarchingCubes(volume, cw, isovalue, true, block, NULL, &vertex_types);
				ITS2TM(oc_its.triangulation, mc_mesh);
			}

			cerr<<"finding connected components"<<endl;
			vector< vector<int> > components;
			GetConnectedComponents(mc_mesh, components);

			vector< vector<Point3> > seed_points;

			for (unsigned i=0; i<components.size(); i++) {

				vector<int> &component = components[i];

				int bound_type = 0;

				for (unsigned v=0; v<component.size(); v++) {
					bound_type |= vertex_types[component[v]];
				}

				// if there are vertices on the boundary of the volume, we'll grow in from them
				if (bound_type & block_type)
					continue;

				// if there are block boundaries on the negative axis, we will grow into this block
				if (bound_type & (MC_VERT_TYPE_XMIN | MC_VERT_TYPE_YMIN | MC_VERT_TYPE_ZMIN))
					continue;

				// we have either a component entirely contained in this block, or it only exits
				// the block on the positive axis sides
				vector<Point3> seeds;
				for (int k=0; k<10; k++) {
					int j = k*component.size() / 10;
					seeds.push_back(mc_mesh.verts[component[j]].point);
				}
				seed_points.push_back(seeds);
			}


			volume.SetBlockSeeds(used_points, used_norms, used_loops, seed_points);
		}




		// project guidance field points
		vector<Point3> gf_pts;
		vector<real_type> gf_ideals;

		vector<Point3> tpts;
		vector<real_type> tideals;

		if (!reuse_guidance) {

			// compute the curvatures of the boundaries for the guidance field
			cerr<<"gf for boundaries"<<endl;
			vector< vector<Point3> > bnd_pts;
			vector< vector<Vector3> > bnd_norms;
			vector<bool> bnd_loops;
			vector< vector<Point3> > cc_seeds;
			volume.GetBlockSeeds(bnd_pts, bnd_norms, bnd_loops, cc_seeds);

			for (int i=0; i<bnd_pts.size(); i++) {
				vector<Point3> &seg = bnd_pts[i];

				for (int j=0; j<seg.size(); j++) {

					int min, max;

					if (bnd_loops[i]) {
						min = j-5;
						max = j+5;
					} else {
						min = std::max(0, j-5);
						max = std::min((int)seg.size()-1, j+5);

						if (min != j-5 || max!=j+5)
							continue;
					}

					vector<Point3> epts;
					vector<real_type> weights;

					for (int k=min; k<=max; k++) {
						epts.push_back(seg[(k+seg.size())%seg.size()]);
						weights.push_back(1);
					}

					real_type k = edge_curvature(epts, weights, j-min);
					gf_pts.push_back(seg[j]);
					gf_ideals.push_back(MaxCurvatureToIdeal(k));
				}
			}


			cerr<<"finding gf points"<<endl;
			ParallelExecutor(idealNumThreads, makeClassFunctor(this,&IsoSurfaceGuidanceField::BuildGuidanceParallel), isovalue, projector, gf_pts, gf_ideals);
			num_orig += gf_pts.size();

			// trim 1
			vector<Point3> tpts1;
			vector<real_type> tideals1;

			vector<int> marked1(gf_ideals.size());
			Trimmer< GetPointVector<Point3> > trimmer1(*this, GetPointVector<Point3>(gf_pts), gf_ideals);
			trimmer1.Trim(marked1);

			for (unsigned i=0; i<marked1.size(); i++) {
				if (!marked1[i]) {
					tpts1.push_back(gf_pts[i]);
					tideals1.push_back(gf_ideals[i]);
				}
			}

			// trim 2
			vector<int> marked(tideals1.size());
			Trimmer< GetPointVector<Point3> > trimmer(*this, GetPointVector<Point3>(tpts1), tideals1);
			trimmer.Trim(marked);

			for (unsigned i=0; i<marked.size(); i++) {
				if (!marked[i]) {
					tpts.push_back(tpts1[i]);
					tideals.push_back(tideals1[i]);
				}
			}


			num_trimmed += tpts.size();
		}

		else {
			cerr<<"using old guidance field!"<<endl;

			extern char *gf_prefix;
			char fname[1024];
			int ab[3];
			Volume::BlockString(gf_prefix==NULL ? "gf" : gf_prefix,
								volume.GetActiveBlock3(ab), fname);
			IsoSurfaceGuidanceField::ReadGuidance(fname, tpts, tideals);

		}

		volume.SetBlockGuidanceField(tpts, tideals);
    }

	volume.SetGuidanceFieldDone();

	if (!reuse_guidance) {
		cerr<<num_orig<<" original guidance field points"<<endl;
		cerr<<num_trimmed<<" guidance field points after trimming"<<endl;
	}
}



IsoSurfaceGuidanceField::~IsoSurfaceGuidanceField() {
}



void* IsoSurfaceGuidanceField::OrderedPointTraverseStart(const Point3 &p) {
    return volume.OrderedPointTraverseStart(p);
}


bool IsoSurfaceGuidanceField::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {
    return volume.OrderedPointTraverseNext(ctx, squared_dist, point, ideal);
}


void IsoSurfaceGuidanceField::OrderedPointTraverseEnd(void *ctx) {
    volume.OrderedPointTraverseEnd(ctx);
}




void IsoSurfaceGuidanceField::Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad) {
#if 0
    pts.resize(kdGetPoint.allpoints.size());
    norms.resize(kdGetPoint.allpoints.size());
    rad.resize(kdGetPoint.allpoints.size());

    for (unsigned i=0; i<kdGetPoint.allpoints.size(); i++) {

		pts[i] = kdGetPoint.allpoints[i];
		if (!volume.NormalAtPoint(pts[i], norms[i])) norms[i] = Vector3(0,0,0);
		rad[i] = ideal_length[i];

    }
#endif
}



void IsoSurfaceGuidanceField::WriteGuidance(const char *fname, const vector<Point3> &pts, const vector<real_type> &ideal) {
#ifdef WIN32
    cerr<<"make sure ./oocgf/ exists!"<<endl;
#else
    mkdir("./oocgf", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    char fullname[1024];
    sprintf(fullname, "./oocgf/%s", fname);
    int size = pts.size();

#ifdef HAS_ZLIB
    gzFile f = gzopen(fullname, "w");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }
		    
    gzwrite(f, &size, sizeof(int));
    gzwrite(f, &pts[0], size*sizeof(Point3));
    gzwrite(f, &ideal[0], size*sizeof(real_type));

    gzclose(f);
#else
    FILE *f = fopen(fullname, "w");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

    fwrite(&size, 1, sizeof(int), f);
    fwrite(&pts[0], size, sizeof(Point3), f);
    fwrite(&ideal[0], size, sizeof(real_type), f);

    fclose(f);
#endif

}


void IsoSurfaceGuidanceField::ReadGuidance(const char *fname, vector<Point3> &pts, vector<real_type> &ideal) {

	//	cerr<<"reading precomputed gf"<<endl;

    char fullname[1024];
    sprintf(fullname, "./oocgf/%s", fname);
    int size;

#ifdef HAS_ZLIB
    gzFile f = gzopen(fullname, "r");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }


    gzread(f, &size, sizeof(int));

    pts.resize(size);
    ideal.resize(size);

    gzread(f, &pts[0], size*sizeof(Point3));
    gzread(f, &ideal[0], size*sizeof(real_type));

    gzclose(f);
#else
    FILE *f = fopen(fullname, "r");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

    fread(&size, 1, sizeof(int), f);

    pts.resize(size);
    ideal.resize(size);

    fread(&pts[0], size, sizeof(Point3), f);
    fread(&ideal[0], size, sizeof(real_type), f);

    fclose(f);
#endif
}


void OOCVolume::WriteSeeds(const char *fname, const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds) {
#ifdef WIN32
    cerr<<"make sure ./oocgf/ exists!"<<endl;
#else
    mkdir("./oocgf", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    char fullname[1024];
    sprintf(fullname, "./oocgf/%s", fname);

#ifdef HAS_ZLIB
    gzFile f = gzopen(fullname, "w");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

	int nboundaries = bnd_pts.size();

    gzwrite(f, &nboundaries, sizeof(int));

	for (int i=0; i<nboundaries; i++) {
		int this_size = bnd_pts[i].size();
		bool this_loop = loops[i];
		gzwrite(f, &this_size, sizeof(int));
		gzwrite(f, &this_loop, sizeof(bool));
		gzwrite(f, &bnd_pts[i][0], sizeof(Point3)*this_size);
		gzwrite(f, &bnd_norms[i][0], sizeof(Vector3)*this_size);
	}


	int nseeds = cc_seeds.size();
    gzwrite(f, &nseeds, sizeof(int));

	for (int i=0; i<nseeds; i++) {
		int this_size = cc_seeds[i].size();
		gzwrite(f, &this_size, sizeof(int));
		gzwrite(f, &cc_seeds[i][0], sizeof(Point3)*this_size);
	}

    gzclose(f);

#else
    FILE *f = fopen(fullname, "w");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

	int nboundaries = bnd_pts.size();

    fwrite(&nboundaries, 1, sizeof(int), f);

	for (int i=0; i<nboundaries; i++) {
		int this_size = bnd_pts[i].size();
		bool this_loop = loops[i];
		fwrite(&this_size, sizeof(int), 1, f);
		fwrite(&this_loop, sizeof(bool), 1, f);
		fwrite(&bnd_pts[i][0], sizeof(Point3)*this_size, 1, f);
		fwrite(&bnd_norms[i][0], sizeof(Vector3)*this_size, 1, f);
	}


	int nseeds = cc_seeds.size();
    fwrite(&nseeds, sizeof(int), 1, f);

	for (int i=0; i<nseeds; i++) {
		int this_size = cc_seeds[i].size();
		fwrite(&this_size, sizeof(int), 1, f);
		fwrite(&cc_seeds[i][0], sizeof(Point3)*this_size, 1, f);
	}

    fclose(f);

#endif

}


void OOCVolume::ReadSeeds(const char *fname, vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) {

	//	cerr<<"reading precomputed gf"<<endl;

    char fullname[1024];
    sprintf(fullname, "./oocgf/%s", fname);


#ifdef HAS_ZLIB
    gzFile f = gzopen(fullname, "r");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

	int nboundaries;
    gzread(f, &nboundaries, sizeof(int));
	bnd_pts.resize(nboundaries);
	bnd_norms.resize(nboundaries);
	loops.resize(nboundaries);

	for (int i=0; i<nboundaries; i++) {
		int this_size;
		bool this_loop;
		gzread(f, &this_size, sizeof(int));
		gzread(f, &this_loop, sizeof(bool));

		bnd_pts[i].resize(this_size);
		gzread(f, &bnd_pts[i][0], sizeof(Point3)*this_size);

		bnd_norms[i].resize(this_size);
		gzread(f, &bnd_norms[i][0], sizeof(Vector3)*this_size);

		loops[i] = this_loop;
	}


	int nseeds;
    gzread(f, &nseeds, sizeof(int));
	cc_seeds.resize(nseeds);
	
	for (int i=0; i<nseeds; i++) {
		int this_size;
		gzread(f, &this_size, sizeof(int));

		cc_seeds[i].resize(this_size);
		gzread(f, &cc_seeds[i][0], sizeof(Point3)*this_size);
	}


    gzclose(f);


#else

    FILE *f = fopen(fullname, "r");
    if (!f) {
		cerr<<"couldn't open "<<fullname<<endl;
		return;
    }

	int nboundaries;
    fread(&nboundaries, sizeof(int), 1, f);
	bnd_pts.resize(nboundaries);
	bnd_norms.resize(nboundaries);
	loops.resize(nboundaries);

	for (int i=0; i<nboundaries; i++) {
		int this_size;
		bool this_loop;
		fread(&this_size, sizeof(int), 1, f);
		fread(&this_loop, sizeof(bool), 1, f);

		bnd_pts[i].resize(this_size);
		fread(&bnd_pts[i][0], sizeof(Point3)*this_size, 1, f);

		bnd_norms[i].resize(this_size);
		fread(&bnd_norms[i][0], sizeof(Vector3)*this_size, 1, f);

		loops[i] = this_loop;
	}


	int nseeds;
    fread(&nseeds, sizeof(int), 1, f);
	cc_seeds.resize(nseeds);
	
	for (int i=0; i<nseeds; i++) {
		int this_size;
		fread(&this_size, sizeof(int), 1,f);

		cc_seeds[i].resize(this_size);
		fread(&cc_seeds[i][0], sizeof(Point3)*this_size, 1, f);
	}


    fclose(f);

#endif
}



//=====================================================================================
//=====================================================================================

Volume::Volume() {
    dim[0]=dim[1]=dim[2]=0;
    aspect[0]=aspect[1]=aspect[2]=1;
    invalidate_all();
    translation = Vector3(0,0,0);
}

Volume::~Volume() {
}

void Volume::SetSplineCoeffs(bool bspline) {
    if (bspline) {
		spline.SetCoefsBSpline();
		cspline.SetCoefsBSpline();
	}
    else {
		spline.SetCoefsCatmullRom();
		cspline.SetCoefsCatmullRom();
	}
}



bool Volume::LocalInfo(const Point3 &p, int cell[3], double xl[3]) const {

    extern bool allow_outside;

    Point3 pt = p - translation;

    for (int i=0; i<3; i++) {
		cell[i] = round_to_minus_inf(pt[i]/aspect[i]);
		xl[i] = pt[i] - aspect[i]*cell[i];
    }

    if (!allow_outside) {
		for (int i=0; i<3; i++) {

			if (cell[i] < 0)
				return false;

			if (cell[i] == dim[i]-1) {
				if (xl[i] > 1e-5)
					return false;
			}

			if (cell[i] > dim[i]-1)
				return false;
		}
    }

    return true;
}



bool Volume::CellCornerValues(const int cell[3], real_type values[2][2][2], real_type isovalue) const {

    double xl[3];
    double nbrs[4][4][4];
    Gather(cell, nbrs);

    double min=nbrs[0][0][0];
    double max=nbrs[0][0][0];
    for (int i=1; i<(4*4*4); i++) {
		min = std::min(min, (&nbrs[0][0][0])[i]);
		max = std::max(max, (&nbrs[0][0][0])[i]);
    }

    // if (min>isovalue || max<isovalue)
	// 	return false;


    for (int x=0; x<2; x++) { 
		xl[0] = x;
		for (int y=0; y<2; y++) {
			xl[1] = y;
			for (int z=0; z<2; z++) {
				xl[2] = z;

				values[x][y][z] = spline.Eval(1, 1, 1, nbrs, xl);
			}
		}
    }

    if (min>isovalue || max<isovalue)
		return false;

    return true;
}


bool Volume::CellSpansValue(const int cell[3], int interior_samples, real_type isovalue) const {

    double xl[3];
    double nbrs[4][4][4];
    Gather(cell, nbrs);

    real_type min = 1e34;
    real_type max =-1e34;

    for (int x=0; x<=interior_samples+1; x++) {
		xl[0] = (real_type)x / (real_type)(interior_samples+1);

		for (int y=0; y<=interior_samples+1; y++) {
			xl[1] = (real_type)y / (real_type)(interior_samples+1);

			for (int z=0; z<=interior_samples+1; z++) {
				xl[2] = (real_type)z / (real_type)(interior_samples+1);

				real_type f = spline.Eval(1, 1, 1, nbrs, xl);

				min = std::min(f,min);
				max = std::max(f,max);

				if (min<=isovalue && max>=isovalue)
					return true;
			}
		}
    }

    return false;
    //    return (min<=isovalue && max>=isovalue);
}


real_type Volume::EvalAtPoint(const Point3 &p) const {

    double xl[3];
    int cell[3];
    double nbrs[4][4][4];

    if (!LocalInfo(p, cell, xl)) return 1e34;
    Gather(cell, nbrs);
    return spline.Eval(GetAspect(0), GetAspect(1), GetAspect(2), nbrs, xl);
}


bool Volume::NormalAtPoint(const Point3 &p, Vector3 &n) const {

    double xl[3];
    int cell[3];
    double nbrs[4][4][4];

    if (!LocalInfo(p, cell, xl)) return false;
    Gather(cell, nbrs);

    double gradient[3];
    spline.Gradient(GetAspect(0), GetAspect(1), GetAspect(2), nbrs, xl, gradient);

    n = Vector3(gradient[0], gradient[1], gradient[2]);
    n.normalize();

    return true;
}


bool Volume::GradientAtPoint(const Point3 &p, Vector3 &g) const {

    double xl[3];
    int cell[3];
    double nbrs[4][4][4];

    if (!LocalInfo(p, cell, xl)) return false;
    Gather(cell, nbrs);

    double gradient[3];
    spline.Gradient(GetAspect(0), GetAspect(1), GetAspect(2), nbrs, xl, gradient);

    g = Vector3(gradient[0], gradient[1], gradient[2]);

    return true;
}



bool Volume::CurvatureAtPoint(const Point3 &p, real_type &k1, real_type &k2) const {

    double xl[3];
    int cell[3];
    double nbrs[4][4][4];

    if (!LocalInfo(p, cell, xl)) return false;
    Gather(cell, nbrs);


    typedef gtb::tmat3<double> mat3;

    mat3 H;
    cspline.Hessian(GetAspect(0), GetAspect(1), GetAspect(2), nbrs, xl, H);

    double gradient[3];
    cspline.Gradient(GetAspect(0), GetAspect(1), GetAspect(2), nbrs, xl, gradient);

    
    gtb::tVector3<double> normal(gradient[0], gradient[1], gradient[2]);
    double gradmag = normal.length();
    normal /= gradmag;

    mat3 nnt(normal[0]*normal, normal[1]*normal, normal[2]*normal);
    mat3 I(gtb::tVector3<double>(1,0,0), gtb::tVector3<double>(0,1,0), gtb::tVector3<double>(0,0,1));
    mat3 P = I - nnt;

    mat3 G = (-1/gradmag) * P*H*P;


    double T = G[0][0] + G[1][1] + G[2][2];	// trace
    double F=0;	// frobenius norm
    for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			F += G[i][j]*G[i][j];
		}
    }
    F = sqrt(F);

    double st = sqrt(std::max(0.0, 2*F*F-T*T));
    k1 = (real_type)((T+st)/2);
    k2 = (real_type)((T-st)/2);

    return true;
}


void Volume::BlockString(const char *prefix, const int b[3], char string[1024]) {

    sprintf(string, "%s.%03d.%03d.%03d.nhdr", prefix, b[0], b[1], b[2]);
}



////////////////////////////////////////////////////////////////////////////////////////////

MultiMaterialVolume::MultiMaterialVolume(const Volume &segmentation) : kdGetPoint(gfpoints)
{

	aspect[0] = segmentation.GetAspect(0);
	aspect[1] = segmentation.GetAspect(1);
	aspect[2] = segmentation.GetAspect(2);

	dim[0] = segmentation.GetDim(0);
	dim[1] = segmentation.GetDim(1);
	dim[2] = segmentation.GetDim(2);

	translation = segmentation.GetTranslation();



    kdtree = NULL;
	cur_mat1=0;
	cur_mat2=1;

	real_type fmin, fmax;
	segmentation.GetBlockMinMax(fmin, fmax);
	if (fmin != (int)fmin ||
		fmax != (int)fmax) {
		cerr<<"range not integers!"<<endl;
	}

	int min = (int)fmin;
	int max = (int)fmax;

	if (min != 0)
		cerr<<"min material id not 0??"<<endl;

	for (int i=min; i<=max; i++) {
		RegularVolume<float> v;
		v.FromSegmentation(segmentation, i);

		RegularVolume<float> *s = new RegularVolume<float>();
		s->Smoothed(v, 1);

		matProbs.push_back(s);
	}
}

MultiMaterialVolume::MultiMaterialVolume(const vector<RegularVolume<float> *> &mats) : kdGetPoint(gfpoints)
{

	aspect[0] = mats[0]->GetAspect(0);
	aspect[1] = mats[0]->GetAspect(1);
	aspect[2] = mats[0]->GetAspect(2);

	dim[0] = mats[0]->GetDim(0);
	dim[1] = mats[0]->GetDim(1);
	dim[2] = mats[0]->GetDim(2);

	translation = mats[0]->GetTranslation();



    kdtree = NULL;
	cur_mat1=0;
	cur_mat2=1;

	matProbs = mats;
}



MultiMaterialVolume::~MultiMaterialVolume() {
	for (int i=0; i<matProbs.size(); i++)
		delete matProbs[i];
    if (kdtree) delete kdtree;
}


real_type MultiMaterialVolume::GetValue(int x, int y, int z) const {
	return (matProbs[cur_mat1]->GetValue(x,y,z) - matProbs[cur_mat2]->GetValue(x,y,z));
}


void MultiMaterialVolume::Gather(const int cell[3], double nbrs[4][4][4]) const {

	double nbrs1[4][4][4];
	double nbrs2[4][4][4];

	matProbs[cur_mat1]->Gather(cell, nbrs1);
	matProbs[cur_mat2]->Gather(cell, nbrs2);

	for (int i=0; i<4*4*4; i++) {
		(&nbrs[0][0][0])[i] = (&nbrs1[0][0][0])[i] - (&nbrs2[0][0][0])[i];
	}
}


int MultiMaterialVolume::NumBlocks() const {
    return 1;
}

void MultiMaterialVolume::SetActiveBlock(int b) const {
}

int MultiMaterialVolume::GetActiveBlock() const {
    return 0;
}


int* MultiMaterialVolume::GetActiveBlock3(int i3[3]) const {
    i3[0] = 0;
    i3[1] = 0;
    i3[2] = 0;
	return i3;
}


void MultiMaterialVolume::GetBlockBox(int min[3], int max[3]) const {
    for (int i=0; i<3; i++) {
		min[i] = 0;
		max[i] = dim[i];
    }
}


void MultiMaterialVolume::GetBlockBox(int block[3], int min[3], int max[3]) const {
    return GetBlockBox(min, max);
}


void MultiMaterialVolume::GetBlockMinMax(real_type &min, real_type &max) const {
	min = -1 * iso_scale;
	max =  1 * iso_scale;
}



void MultiMaterialVolume::WriteValue(int x, int y, int z, FILE *f) const {
	float v = GetValue(x,y,z);
    fwrite(&v, 1, sizeof(float), f);
}


#ifdef HAS_ZLIB
void MultiMaterialVolume::WriteValue(int x, int y, int z, gzFile f) const {
	float v = GetValue(x,y,z);
    gzwrite(f, &v, sizeof(float));
}
#endif

const char* MultiMaterialVolume::TypeString() const {
    return "float";
}


void MultiMaterialVolume::SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal) {

	if (kdtree) 
		delete kdtree;
			
    gfpoints = pts;
    ideal_lengths = ideal;

	// cant make empty trees!
	if (pts.size() > 0) {
		kdtree = new kdtree_type(10, Box3::bounding_box(pts), kdGetPoint);
		for (unsigned i=0; i<gfpoints.size(); i++)
			kdtree->Insert(i);
		kdtree->MakeTree();
	} else {
		kdtree = NULL;
	}
}


void MultiMaterialVolume::SetGuidanceFieldDone() {
}


void MultiMaterialVolume::SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds) {

    // boundary_pts = bnd_pts;
    // boundary_norms = bnd_norms;
    // boundary_loops = loops;
    // seed_pts = cc_seeds;
}

void MultiMaterialVolume::GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const {

    // boundary_pts = bnd_pts;
    // boundary_norms = bnd_norms;
    // boundary_loops = loops;
    // seed_pts = cc_seeds;
}


void MultiMaterialVolume::GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {

    int bmins[3], bmaxs[3];
    GetBlockBox(bmins, bmaxs);

    box = Box3(Point3(1e34,1e34,1e34), Point3(-1e34,-1e34,-1e34));
    for (int i=0; i<3; i++) {
		if (GetAspect(i) > 0) {
			box.set_min(i, bmins[i]*GetAspect(i) + translation[i]);
			box.set_max(i, bmaxs[i]*GetAspect(i) + translation[i]);
		} else {
			box.set_max(i, bmins[i]*GetAspect(i) + translation[i]);
			box.set_min(i, bmaxs[i]*GetAspect(i) + translation[i]);
		}
    }

    // pts = boundary_pts;
    // norms = boundary_norms;
    // loops = boundary_loops;
    // seeds = seed_pts;
}


bool MultiMaterialVolume::GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const {
    return false;
}


int MultiMaterialVolume::GetBlockForPoint(const Point3 &p) const {
	return 0;
}


void* MultiMaterialVolume::OrderedPointTraverseStart(const Point3 &p) {
    if (kdtree) {
		return new kdtree_type::OrderedIncrementalTraverse(*kdtree, p);
	} else {
		return NULL;
	}
}



bool MultiMaterialVolume::OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) {

    if (!ctx) {
		return false;
	}

	kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
	if (kdOrderedTraverse->empty()) return false;
    int n = kdOrderedTraverse->GetNext(squared_dist);
    point = gfpoints[n];
    ideal = ideal_lengths[n];
    return true;
}



void MultiMaterialVolume::OrderedPointTraverseEnd(void *ctx) {
	if (ctx) {
		kdtree_type::OrderedIncrementalTraverse *kdOrderedTraverse = (kdtree_type::OrderedIncrementalTraverse*)ctx;
		delete kdOrderedTraverse;
	}
}



int MultiMaterialVolume::MaterialAtPoint(const Point3 &p) const {

	real_type bestv = matProbs[0]->EvalAtPoint(p);
	int besti = 0;

	for (int i=1; i<matProbs.size(); i++) {
		real_type tv = matProbs[i]->EvalAtPoint(p);
		if (tv > bestv) {
			bestv = tv;
			besti = i;
		}
	}

	return besti;
}



void MultiMaterialVolume::SetSplineCoeffs(bool bspline) {
	for (int i=0; i<matProbs.size(); i++) {
		matProbs[i]->SetSplineCoeffs(bspline);
	}

	if (bspline)
		spline.SetCoefsBSpline();
	else
		spline.SetCoefsCatmullRom();

	// always use bspline for curvature here!
	cspline.SetCoefsBSpline();
}
