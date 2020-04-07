
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



#ifndef _TRIANGULATE_ISO
#define	_TRIANGULATE_ISO

#include "interval.h"
#include "spline.h"
#include "sfcurve.h"


template <typename T>
	class Array3 {
	public:
	Array3() { }
	Array3(const T &a, const T &b, const T &c) {
		x[0] = a;
		x[1] = b;
		x[2] = c;
	}
	T& operator[](int i) { return x[i]; }
	const T& operator[](int i) const { return x[i]; }
	bool operator==(const Array3 &rhs) const {
		return (x[0]==rhs[0] && 
				x[1]==rhs[1] && 
				x[2]==rhs[2]);
	}
	bool operator<(const Array3 &rhs) const {
		for (int i=0; i<3; i++) {
			if (x[i] < rhs[i])
				return true;
			if (x[i] > rhs[i])
				return false;
		}
		return false;
	}
	T x[3];
};



// abstraction of volume interface we need for standard volumes and for out of core volumes
class Volume : public gtb::tModel<real_type> {
    public:

    Volume();
    virtual ~Volume();

    // operations we can perform
    virtual void SetSplineCoeffs(bool bspline);
    bool CellSpansValue(const int cell[3], int interior_samples, real_type isovalue) const;
    bool CellCornerValues(const int cell[3], real_type values[2][2][2], real_type isovalue) const;
    real_type EvalAtPoint(const Point3 &p) const;
    bool NormalAtPoint(const Point3 &p, Vector3 &n) const;
    bool GradientAtPoint(const Point3 &p, Vector3 &n) const;
    bool CurvatureAtPoint(const Point3 &p, real_type &k1, real_type &k2) const;
    bool LocalInfo(const Point3 &p, int cell[3], double xl[3]) const;
    int GetDim(int i) const { return dim[i]; }
    real_type GetAspect(int i) const { return aspect[i]; }
    Vector3 GetTranslation() const { return translation; }



    // abstract interface
    virtual void Gather(const int cell[3], double nbrs[4][4][4]) const = 0;
    virtual real_type GetValue(int x, int y, int z) const = 0;

    virtual int NumBlocks() const = 0;

    // THERE MUST NOT BE ANY BLOCK ACCESSES CONCURRENT WITH SetActiveBlock()!!!!
    virtual void SetActiveBlock(int b) const = 0;
    virtual int GetActiveBlock() const = 0;
	virtual int* GetActiveBlock3(int i3[3]) const = 0;

    virtual void GetBlockBox(int block[3], int min[3], int max[3]) const = 0;
    virtual void GetBlockBox(int min[3], int max[3]) const = 0;
    virtual void GetBlockMinMax(real_type &min, real_type &max) const = 0;

    static void BlockString(const char *prefix, const int b[3], char string[1024]);

    virtual void WriteValue(int x, int y, int z, FILE *f) const = 0;
#ifdef HAS_ZLIB
    virtual void WriteValue(int x, int y, int z, gzFile f) const = 0;
#endif
    virtual const char* TypeString() const = 0;


    virtual void SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal) = 0;
    virtual void SetGuidanceFieldDone() = 0;

    virtual void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3 > > &seeds) const = 0; // THIS CALLS SetActiveBlock()!!!
    virtual bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3 > > &seeds) const = 0; // THIS CALLS SetActiveBlock()!!!
	virtual int  GetBlockForPoint(const Point3 &p) const = 0;
    virtual void SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds) = 0;
    virtual void GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const = 0;
	virtual int GetNumBlockLoads() const = 0;

    virtual void* OrderedPointTraverseStart(const Point3 &p) = 0;
    virtual bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal) = 0;
    virtual void OrderedPointTraverseEnd(void *ctx) = 0;



    // Model stuff
    void compute_bounding_box() const {
		Point3 minp = Point3(0,0,0);
		Point3 maxp = Point3(GetAspect(0)*(GetDim(0)-1),
							 GetAspect(1)*(GetDim(1)-1),
							 GetAspect(2)*(GetDim(2)-1));
		minp += GetTranslation();
		maxp += GetTranslation();

		_bounding_box = Box3(Point3(std::min(minp[0], maxp[0]),
									std::min(minp[1], maxp[1]),
									std::min(minp[2], maxp[2])),
							 Point3(std::max(minp[0], maxp[0]),
									std::max(minp[1], maxp[1]),
									std::max(minp[2], maxp[2])));

		_is_bounding_box_valid = true;
    }


    void compute_centroid() const {
		_centroid = bounding_box().centroid();
		_is_centroid_valid = true;
    }

    TrivariateSpline<double> spline;
    TrivariateSpline<double> cspline; // for curvature

    protected:
    int dim[3];
    real_type aspect[3];
    Vector3 translation;

};



class NHDR {
    public:
    enum SourceType {
		Float,
		UChar,
		Short,
		Unknown,
    };

    SourceType type;
    int sizes[3];
    real_type spacings[3];
    bool endian;
    Vector3 translation;
    char datafile[1024];
};


// a function sampled on a regular grid
template <typename T>
	class RegularVolume : public Volume {

    public:
    RegularVolume();
    ~RegularVolume();

    bool ReadNRRD(const NHDR &nhdr);
    void FromSegmentation(const Volume &v, int mat);
	void FromMultiMaterial(vector<Volume*> &vols, int i);
    void Smoothed(const Volume &v, int width);

    void Gather(const int cell[3], double nbrs[4][4][4]) const;
    real_type GetValue(int x, int y, int z) const;
    int NumBlocks() const;
    void SetActiveBlock(int b) const;
    int GetActiveBlock() const;
	int* GetActiveBlock3(int i3[3]) const;
    void GetBlockBox(int block[3], int min[3], int max[3]) const;
    void GetBlockBox(int min[3], int max[3]) const;
    void GetBlockMinMax(real_type &min, real_type &max) const;

    void WriteValue(int x, int y, int z, FILE *f) const;
#ifdef HAS_ZLIB
    void WriteValue(int x, int y, int z, gzFile f) const;
#endif
    const char* TypeString() const;

    void SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal);
    void SetGuidanceFieldDone();

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
	int  GetBlockForPoint(const Point3 &p) const;
    void SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds);
    void GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const;
	int GetNumBlockLoads() const { return 0; }

    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);


    private:
    bool ReadSource(const char *fname, bool bigendian=false);

    T *data;
    T min, max;

    // guidance field stuff
    GetPointVector<Point3> kdGetPoint;
    vector<Point3> gfpoints;
    vector<real_type> ideal_lengths;
    typedef gtb::KDTree<int, real_type, GetPointVector<Point3> > kdtree_type;
    kdtree_type *kdtree;

    // stored boundary info
    vector< vector<Point3> > boundary_pts;
    vector< vector<Vector3> > boundary_norms;
    vector<bool> boundary_loops;
    vector< vector<Point3> > seed_pts;
};


class MultiMaterialVolume : public Volume {

    public:
    MultiMaterialVolume(const Volume &segmentation);
    MultiMaterialVolume(const vector<RegularVolume<float> *> &mats);
    ~MultiMaterialVolume();


    virtual void SetSplineCoeffs(bool bspline);

    void Gather(const int cell[3], double nbrs[4][4][4]) const;
    real_type GetValue(int x, int y, int z) const;
    int NumBlocks() const;
    void SetActiveBlock(int b) const;
    int GetActiveBlock() const;
	int* GetActiveBlock3(int i3[3]) const;
    void GetBlockBox(int block[3], int min[3], int max[3]) const;
    void GetBlockBox(int min[3], int max[3]) const;
    void GetBlockMinMax(real_type &min, real_type &max) const;

    void WriteValue(int x, int y, int z, FILE *f) const;
#ifdef HAS_ZLIB
    void WriteValue(int x, int y, int z, gzFile f) const;
#endif
    const char* TypeString() const;

    void SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal);
    void SetGuidanceFieldDone();

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
	int  GetBlockForPoint(const Point3 &p) const;
    void SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds);
    void GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const;
	int GetNumBlockLoads() const { return 0; }

    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);


	int MaterialAtPoint(const Point3 &p) const;

	int cur_mat1, cur_mat2;
	vector< RegularVolume<float>* > matProbs;

    private:

    // guidance field stuff
    GetPointVector<Point3> kdGetPoint;
    vector<Point3> gfpoints;
    vector<real_type> ideal_lengths;
    typedef gtb::KDTree<int, real_type, GetPointVector<Point3> > kdtree_type;
    kdtree_type *kdtree;
};



// out of core regular grid volume
class OOCVolume : public Volume {

    public:
    OOCVolume();
    ~OOCVolume();

    bool Read(const char *fname);

    void Gather(const int cell[3], double nbrs[4][4][4]) const;
    real_type GetValue(int x, int y, int z) const;
    int NumBlocks() const;
    void SetActiveBlock(int b) const;
    int GetActiveBlock() const;
	int* GetActiveBlock3(int i3[3]) const;
    void GetBlockBox(int block[3], int min[3], int max[3]) const;
    void GetBlockBox(int min[3], int max[3]) const;
    void GetBlockMinMax(real_type &min, real_type &max) const;


    void WriteValue(int x, int y, int z, FILE *f) const;
#ifdef HAS_ZLIB
    void WriteValue(int x, int y, int z, gzFile f) const;
#endif
    const char* TypeString() const;

    void SetBlockGuidanceField(const vector<Point3> &pts, const vector<real_type> &ideal);
    void SetGuidanceFieldDone();

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3> > &seeds) const;
	int  GetBlockForPoint(const Point3 &p) const;
    void SetBlockSeeds(const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds);
    void GetBlockSeeds(vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds) const;
	int GetNumBlockLoads() const { return num_block_loads; }

    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);



    private:

    static void WriteSeeds(const char *fname, const vector< vector<Point3> > &bnd_pts, const vector< vector<Vector3> > &bnd_norms, const vector<bool> &loops, const vector< vector<Point3> > &cc_seeds);
    static void ReadSeeds(const char *fname, vector< vector<Point3> > &bnd_pts, vector< vector<Vector3> > &bnd_norms, vector<bool> &loops, vector< vector<Point3> > &cc_seeds);


    char block_prefix[1024];
    class BlockData {
		public:
		BlockData() : rv(NULL), load_stamp(-1), ordered_traverse_refcount(0) { }
		~BlockData() { if (rv) delete rv; }

		Volume *rv; // really RegularVolume<>'s
		int64 load_stamp;
		int ordered_traverse_refcount;
    };


    int blocksize[3];
    int nblocks[3];

    // must be thread safe!
    void BlockString(const int b[3], char string[1024]) const;
    void BlockDim(const int b[3], int d[3]) const;

    bool BlockInActiveRegion(const int block[3]) const;
    BlockData& BeginBlockAccess(const int block[3]) const;
    void EndBlockAccess(const int block[3]) const;


    mutable thlib::CSObject cs;
    mutable vector< vector< vector<BlockData> > > blocks;
    mutable volatile int64 blockToLoadASAP;

	int num_block_loads;
	

    mutable int active_index;
    mutable int active_block[3];
    mutable int64 active_zindex;

    mutable volatile bool active_region_loading;

	// we use a condition variable to wake up the fetching thread
	mutable volatile int need_fetch_work; // 0=no work, 1=work, 2=kill
#ifndef WIN32
	mutable pthread_mutex_t fetch_mutex;
	mutable pthread_cond_t fetch_condition;
#endif
	bool FetcherWait();
	void WakeupFetcher(bool kill=false) const;


    void DrawBlock(const Array3<int> &b, const Point3 &color) const;
    void DrawLoaded(std::set<Array3<int> > &prefetched,
					std::set<Array3<int> > &grab_bag,
					vector< Array3<int> > &prefetch_list) const;
    void ValidateFetchState(std::set<Array3<int> > &prefetched,
							std::set<Array3<int> > &grab_bag,
							vector< Array3<int> > &prefetch_list) const;

    thlib::Thread *prefetcher;
    int64 load_counter;

    void* PrefetchMain();
    static void* PrefetchMain_s(void *t) {
		return ((OOCVolume*)t)->PrefetchMain();
    }
    void LoadBlock(const int b[3]);
    void RemoveLRU(std::set< Array3<int> > &bag);
  

    class OrderedTraverseInfo {
		public:

		bool operator>(const OrderedTraverseInfo &rhs) const {
			return (neg_dist_squared>rhs.neg_dist_squared);
		}
		bool operator<(const OrderedTraverseInfo &rhs) const {
			return (neg_dist_squared<rhs.neg_dist_squared);
		}

		real_type neg_dist_squared;
		int ring;

		Array3<int> block;
		Point3 p;
		real_type ideal;
    };

	class OrderedTraverseContext {
		public:
		gtb::fast_pq<OrderedTraverseInfo>  pq;
		std::map< Array3<int>, void*> subctx;
		Point3 p;
	};

    void GetShellBlocks(int block[3], int shell, vector< Array3<int> > &shell_blocks) const;
    real_type DistToRingBoundary(const Point3 &p, int ring) const;
    real_type DistToBlock(int block[3], const Point3 &p) const;

    volatile bool guidance_computed;
};

Volume* VolumeRead(const char *fname);




class IsoSurfaceProjector : public SurfaceProjector {
    public:
    IsoSurfaceProjector(const Volume &vol, real_type iso);
    virtual ~IsoSurfaceProjector();
    int ProjectPoint(const FrontElement &base1, const FrontElement &base2, const Point3 &fp, const Vector3 &fn, Point3 &tp, Vector3 &tn, real_type *curvature) const;
    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn) const;
    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn, const Vector3 *plane) const;
    int ProjectPoint(const Point3 &fp, Point3 &tp, Vector3 &tn, const Point3 &center, const Vector3 &axis) const;

    void GetFirstBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3 > > &seeds) const {
		volume.GetFirstBlock(box, pts, norms, loops, seeds);
    }
    bool GetNextBlock(Box3 &box, vector< vector<Point3> > &pts, vector< vector<Vector3> > &norms, vector<bool> &loops, vector< vector<Point3 > > &seeds) const {
		return volume.GetNextBlock(box, pts, norms, loops, seeds);
    }
	int  GetBlockForPoint(const Point3 &p) const {
		return volume.GetBlockForPoint(p);
	}


    real_type GetIsoValue() const { return isovalue; }

    private:

    real_type isovalue;
    const Volume &volume;



    class RingEval {

		public:
		RingEval(const IsoSurfaceProjector &t, const Point3 &c, const Vector3 &u, const Vector3 &v) :
			isp(t), center(c), udir(u), vdir(v) { last_cell[0] = last_cell[1] = last_cell[2] = -1; }

		Point3 point(real_type theta) const {
			return center + (real_type)cos(theta)*udir + (real_type)sin(theta)*vdir;
		}

		real_type operator()(real_type theta) {
			Point3 x = point(theta);
	    
			int this_cell[3];
			double xl[3];	// local coordinate

			isp.volume.LocalInfo(x, this_cell, xl);
	    
			if (this_cell[0]!=last_cell[0] || this_cell[1]!=last_cell[1] || this_cell[2]!=last_cell[2]) {
	      
				isp.volume.Gather(this_cell, nbrs);
	      
				last_cell[0]=this_cell[0];
				last_cell[1]=this_cell[1];
				last_cell[2]=this_cell[2];
			}


			double e = isp.volume.spline.Eval(isp.volume.GetAspect(0),
											  isp.volume.GetAspect(1),
											  isp.volume.GetAspect(2),
											  nbrs, xl);

			return (real_type)(e - isp.isovalue);
		}

        
		const IsoSurfaceProjector &isp;
		int last_cell[3]; 
		double nbrs[4][4][4];

		Point3 center;
		Vector3 udir, vdir;
    };
};






class IsoSurfaceGuidanceField : public GuidanceField {
    public:

    IsoSurfaceGuidanceField(IsoSurfaceProjector &projector, Volume &vol,
							real_type rho, real_type min_step, real_type max_step, real_type reduction, int block=-1);
    IsoSurfaceGuidanceField(IsoSurfaceProjector &projector, MultiMaterialVolume &vol, vector< vector<Point3> > &loops, const vector< vector<int> > &loop_mats,
							real_type rho, real_type min_step, real_type max_step, real_type reduction);
    ~IsoSurfaceGuidanceField();

    // we basically implement all this in the Volume class, so this stuff just passes it along
    void* OrderedPointTraverseStart(const Point3 &p);
    bool OrderedPointTraverseNext(void *ctx, real_type &squared_dist, Point3 &point, real_type &ideal);
    void OrderedPointTraverseEnd(void *ctx);

    void Extract(vector<Point3> &pts, vector<Vector3> &norms, vector<real_type> &rad);

    static void WriteGuidance(const char *fname, const vector<Point3> &pts, const vector<real_type> &ideal);
    static void ReadGuidance(const char *fname, vector<Point3> &pts, vector<real_type> &ideal);


    private:

    void BuildGuidanceParallel(int nt, int id, CSObject &cs, real_type isovalue, IsoSurfaceProjector &projector, vector<Point3> &pts, vector<real_type> &ideals) const;
    void BuildMultiGuidanceParallel(int nt, int id, CSObject &cs, MultiMaterialVolume &volume, IsoSurfaceProjector &projector, vector<Point3> &pts, vector<real_type> &ideals) const;


	Interval BoundSpline(const double nbrs[4][4][4], const double coefs[4][4][4], const int partial[3], const Interval xspan_s[3], int dx_sign, int dy_sign, int dz_sign) const;
	void SampleCell(int id, int cell[3], double nbrs[4][4][4], double coefs[4][4][4],
					const Interval xspan_s[3], IsoSurfaceProjector &projector, 
					vector<Point3> &pts, vector<real_type> &ideals, 
					const int third_sign[3][3][3],
					int depth=0) const;

	struct a444_wrap { real_type n[4][4][4]; };
	void SampleCell(int id, int cell[3], vector<a444_wrap> &nbrs, vector<a444_wrap> &coefs, int considerbits,
					const Interval xspan_s[3], IsoSurfaceProjector &projector, MultiMaterialVolume &volume,
					vector<Point3> &pts, vector<real_type> &ideals, 
					int depth=0) const;
	

	int BoundaryRecursionDepth(int axis, double nbrs[4][4][4], double coefs[4][4][4], real_type isovalue, const Interval xspan[3], int depth) const;
	void GetBlockFaceBoundaries(real_type isovalue, int axis, bool min, vector< vector<Point3> > &bsegments, vector<bool> &bloops) const;
	bool AtMostSingleIntersection(double nbrs[4][4][4], const Point3 &x1, const Point3 &x2, real_type isovalue) const;


    Volume &volume;
};




#define MC_VERT_TYPE_XMIN 0x01
#define MC_VERT_TYPE_XMAX 0x02
#define MC_VERT_TYPE_YMIN 0x04
#define MC_VERT_TYPE_YMAX 0x08
#define MC_VERT_TYPE_ZMIN 0x10
#define MC_VERT_TYPE_ZMAX 0x20

void MarchingCubes(const Volume &v, TriangulatorController &tc, real_type isovalue, bool rootfinding=false, int mc_block=-1, vector<int64> *intersected_cells=NULL, vector<int> *vertex_types=NULL);

void MarchingCubesStats(const Volume &v, TriangulatorController *tc, real_type isovalue,
						real_type &cellcount, real_type &wcellcount,
						real_type &tricount, real_type &wtricount,
						real_type &area, real_type &warea);



#endif
