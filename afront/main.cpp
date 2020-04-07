
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
#include <stdlib.h>
#include <ThreadLib/threadslib.h>
#include <gtb/ui/commandline.h>
#include <viewer/viewer.h>
#include "uheap.h"
#include "front.h"
#include "parallel.h"
#include "guidance.h"
#include "triangulator.h"
#include "triangulate_mesh.h"
#include "triangulate_mls.h"
#include "triangulate_iso.h"
#include "triangulate_tet.h"
#include "triangulate_csg.h"

#include "output_controller_obj.h"
#include "output_controller_hhm.h"
#include "output_controller_sma.h"
#include "output_controller_gui.h"
#include "edgeflipper.h"
#include "crease.h"

#include "FLF_io.h"
#include "PC_io.h"

#include "interval.h"



using namespace std;

#ifdef __x86_64__
// this used to be 4, but it appears to be about 10-15% faster with 8
// would we be better off with 2 theads on a uniprocessor?
int idealNumThreads = 8;
#else
int idealNumThreads = 1;
#endif



OutputControllerGui *gui=NULL;
static char command_line[2048];


// can be modified by the worker thread, but only inside of the draw mutex.	 Draw thread will never modify it
static TriangleMesh meshes[2];
static Volume *volume = NULL;
static surfel_set points;
static TetMesh tetmesh;
static TriangleMesh tetshell;

// const versions for the drawing thread so it can't accidently screw it up
const TriangleMesh *cmeshes[2] = { &meshes[0], &meshes[1] };
const Volume *cvolume = NULL;
const surfel_set *cpoints = &points;
const TetMesh *ctetmesh = &tetmesh;
const TriangleMesh *ctetshell = &tetshell;

char defoutname[] = "outmesh.m";
char *outname = defoutname;
char *gf_prefix = NULL;
char *histname = NULL;
char *histfigname = NULL;
static char *gffile = NULL;
real_type rho = (real_type)(M_PI * 2.0 / 32.0); // ~ 0.2
static real_type sharp = -2.0; // if cos dihedralangle < sharp, treat it as crease
static real_type min_step = (real_type) 0;
static real_type max_step = (real_type) 1000000;
real_type reduction = (real_type)0.8;
real_type radius_factor = 2.0f;
static bool csg_guidance_blend = true;
static bool failsafe = true;	
static int small_crease = 20;
static int adamson = 0;
static int nielson = 0;
static bool draw_messages = true;
int boundary_width = 2;
real_type fence_scale = 1.0;
int curvature_sub = 1;
real_type curvature_spread = 2;
int eval_sub = 3;
bool trim_guidance = true;
int mls_gf_fix = 0;
real_type value_tol = 1e-8;
real_type step_tol = 1e-4;
int max_iter = 50;
int trim_bin_size = 1000000;
real_type bad_connect_priority = 1e34;
bool grow_from_connect = true;
real_type boundary_dist = 0;
real_type noise_threshold = 0;
bool rmls = false;
int rmls_knn = 120;
bool stream_mc = false;
char *meshfrontfile = NULL;

int max_grab_bag = 5;
int prefetch_rings = 1;
int prefetch_lookahead = 1;

bool feature_interference_hack = false;
extern bool reuse_guidance;
extern bool reuse_seeds;
bool retro = false;

extern real_type iso_scale;
static int mc_stats=0;

bool is_segmentation = false;

GuidanceField *guidance = NULL;
static ControllerWrapper *controller = NULL;
static Triangulator *triangulator = NULL;

thlib::CSObject *critical_section;


bool noedgeflip = false;
bool endswith(const char *s, const char *end) {
	int slen = strlen(s);
	int endlen = strlen(end);
	if (slen<endlen) return false;
	return (!stricmp(&s[slen-endlen], end));
}


void fileappend(const char *fname, float num) {
	char fullname[1024];
	FILE *f = fopen(strcat(strcpy(fullname, outname),fname) , "ab");
	fprintf(f, "%g\n", num);
	fclose(f);
}


// keep these here for drawing them if we want
vector<OutputControllerEdgeFlipper*> edgeflippers;
OutputController* GetOutputController(bool edgeflip=true, bool ocgui=true) {

	OutputController *ret = NULL;

	// add the edge flipper at the beginning if we want it
	if (edgeflip && !noedgeflip) {
		edgeflippers.clear();

		edgeflippers.push_back(new OutputControllerEdgeFlipper());
		OutputController::AddControllerToBack(ret, edgeflippers.back());

		edgeflippers.push_back(new OutputControllerEdgeFlipper());
		OutputController::AddControllerToBack(ret, edgeflippers.back());
	}

	if (endswith(outname, ".m")) {
		OutputController::AddControllerToBack(ret, new OutputControllerHHM(outname));
	} else if (endswith(outname, ".obj")) {
		OutputController::AddControllerToBack(ret, new OutputControllerOBJ(outname));
	} else if (endswith(outname, ".sma") ||
			   endswith(outname, ".smb") ||
			   endswith(outname, ".smc") ||
			   endswith(outname, ".smd") ||
			   endswith(outname, ".sma.gz") ||
			   endswith(outname, ".smb.gz") ||
			   endswith(outname, ".smc.gz") ||
			   endswith(outname, ".smd.gz")) {
		OutputController::AddControllerToBack(ret, new OutputControllerSMA(outname));
	} else {
		cerr<<"unknow output file type: "<<outname<<endl;
	}

	if (ocgui && gui) {
		OutputController::AddControllerToBack(ret, gui);
	}

	return ret;
}


bool allow_outside=false;


void redrawAll(){
	if (gui) {
		gui->Extract(triangulator);
		gui->Redraw();
	}
}

void redrawAndWait(int key, bool force){
	if(!gui) return;
	redrawAll();

	if (force || gui->GetWaitForUser()) {

		/*
		dbgClear();
		for (int i=0; i<edgeflippers.size(); i++) {
			edgeflippers[i]->Draw();
		}
		*/


		char message[1024];
		sprintf(message, "waiting for key <%c>\n", key);
		if (draw_messages)
			gui->SetMessage(0, message);
		redrawAll();
		gui->WaitForKey(key);
		gui->SetMessage(0, NULL);
	}
}

void waitForKey(int key) {
	if(!gui) return;
	if (gui->GetWaitForUser())
		gui->WaitForKey(key);
}


int GetUILastKey() {
	if (!gui) return 0;
	return gui->LastKey();
}

bool getUISelectPoint(Point3 &p, int &win) {

	if (!gui) return false;

	int button;
	gui->LastClick(p[0], p[1], p[2], win, button);
	return (win>=0 && button==GLUT_RIGHT_BUTTON);
}




void PerpVectors(const Vector3 &n, Vector3 &udir, Vector3 &vdir) {

	int m = (fabs(n[0])<fabs(n[1])) ? 0 : 1;
	if (fabs(n[2])<fabs(n[m])) m=2;

	udir=Vector3(0,0,0);
	udir[m] = 1;

	vdir = udir.cross(n); vdir.normalize();
	udir = vdir.cross(n); udir.normalize();
}



/////////////////////////////////////////////////////////////////////////////////

void do_gui() {

	if (gui) {
		cerr<<"gui already started!"<<endl;
		return;
	}

	gui = new OutputControllerGui(command_line);
	gui->Start();

	critical_section = &gui->GetCriticalSection();

	redrawAll();
	gui->ViewAll();
	gui->SetUnsharpWidth(0);
}


void do_quit(){
	exit(0);
}

//////////////////////////////////////////////
// the main functions that get things going

void extract_boundary_loops(const TriangleMesh &m, vector< vector<int> > &loops) {
	loops.resize(0);

	vector<bool> needadd(m.verts.size());
	for (unsigned i=0; i<needadd.size(); i++) 
		needadd[i] = false;

	for (unsigned v=0; v<m.verts.size(); v++) {
		if (m.verts[v].someface >= 0) {
			TriangleMesh::VertexVertexIteratorI vvi(m, v);
			needadd[v] = vvi.isonedge();
		}
	}


	while (1) {

		int start = -1;
		for (unsigned i=0; i<needadd.size(); i++) {
			if (needadd[i]) {
				start = i;
				break;
			}
		}

		if (start == -1) break;


		vector<int> tloop;
		tloop.push_back(start);
		needadd[start] = false;

		while (1) {
			TriangleMesh::VertexVertexIteratorI vvi(m, tloop.back());
			int next = *vvi;
			if (next == start) break;
			needadd[next] = false;
			tloop.push_back(next);
		}
		
		loops.push_back(tloop);
	}
}


void compute_pointset_normals_parallel(int nt, int id, CSObject &cs, surfel_set &ss, CProjection &proj) {
	
	for (unsigned i=id; i<ss.size(); i+=nt) {
		Point3 r1;
		Vector3 n;
		proj.PowellProject(ss.vertex(i), r1, n);
		ss.normal(i) = n;
	}
}

void orient_normals(int startv, surfel_set &ss, CProjection &proj, bool all) {
	// orient the normals
	vector<bool> set(ss.size());
	for (unsigned i=0; i<ss.size(); i++) set[i]=false;

	
	gtb::fast_pq< std::pair<real_type, std::pair<int,int> > > pq;
	
	pq.push(std::pair<real_type, std::pair<int,int> > (0, std::pair<int,int>(-1,startv) ));

	bool first=true;
	int numflipped = 0;

	while (!pq.empty()) {

		std::pair<real_type, std::pair<int,int> > top = pq.top();
		pq.pop();

		if (!set[top.second.second]) {

			if (first) {
				set[top.second.second] = true;
			}

			surfelset_view nbhd(ss);
			proj.extract(ss.vertex(top.second.second), ((real_type)1)*proj.point_radius(ss.vertex(top.second.second)), nbhd);

			Vector3 setnorm(0,0,0);
			for (unsigned i=0; i<nbhd.size(); i++) {
				if (!set[nbhd.get_index(i)]) continue;
				setnorm += ss.normal(nbhd.get_index(i));
			}
			setnorm.normalize();

			if (setnorm.dot(ss.normal(top.second.second)) < 0) {
				ss.normal(top.second.second).flip(); 
				numflipped++;
				if (!all && numflipped > 1000)
					break;
			} else {
				if (!all) {
					if (first) {
						first = false;
					} else {
						set[top.second.second] = true;
						continue;
					}					
				}
			}

			set[top.second.second] = true;

			for (unsigned i=0; i<nbhd.size(); i++) {
				if (set[nbhd.get_index(i)]) continue;

				if (all)
					pq.push(std::pair<real_type, std::pair<int,int> > (top.first - Point3::distance(ss.vertex(top.second.second),nbhd.vertex(i)), std::pair<int,int>(top.second.second,nbhd.get_index(i)) ));
				else {
					if (ss.normal(startv).dot(ss.normal(nbhd.get_index(i)))<-0.5)
						pq.push(std::pair<real_type, std::pair<int,int> > (top.first - Point3::distance(ss.vertex(top.second.second),nbhd.vertex(i)), std::pair<int,int>(top.second.second,nbhd.get_index(i)) ));
					else
						set[nbhd.get_index(i)] = true;
				}
			}
		}
	}
}


void compute_pointset_normals(surfel_set &ss, CProjection &proj) {

	ss.normals().resize(ss.size());
	ParallelExecutor(idealNumThreads, compute_pointset_normals_parallel, ss, proj);
	orient_normals(0, ss, proj, true);
}

void do_compute_pointset_normals() {
	SmoothMLSProjector psp(points, adamson);
	psp.SetRadiusFactor(radius_factor);
	compute_pointset_normals(points, psp._projector);
}


int do_tri_mesh(int argc, char* argv[]) {
	
	int ret=1;
	int csubdiv=0;
	//	char *command = argv[1];
	if (argc>1 && argv[1][0]!='-') {
		csubdiv = atoi(argv[1]);
		ret++;
	}



	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;

	const TriangleMesh &mesh = *cmeshes[0];
	MeshProjector projector(mesh);
	guidance = new MeshGuidanceField(csubdiv, mesh, rho, min_step, max_step, reduction);
	controller = new ControllerWrapper(guidance, &projector, GetOutputController());


	vector< vector<int> > boundaries;
	vector< Point3 > crease_points;
	vector< Vector3 > crease_onormals;
	vector< vector<int> > crease_indices;
	vector< vector<Vector3> > crease_normals;
	vector< int > corner_triangles;

	vector< vector<Point3> > ipts;
	vector< vector<Vector3> > inorms;


#if 1  // this kills pensatore for some reason
	
	if (boundary_dist != 0) {
		projector.FindSoupBoundaries(boundaries, boundary_dist);
	} else if (meshfrontfile) {
		FILE *f = fopen(meshfrontfile, "r");
		if (!f) {
			cerr<<"couldn't open mesh front file "<<meshfrontfile<<endl;
		} else {

			vector<Point3> ploop;
			vector<Vector3> nloop;

			while (!feof(f)) {
				char line[1000];
				if (fgets(line, sizeof(line), f) == 0) break;

				double x,y,z;
				sscanf(line,"%lf %lf %lf", &x, &y, &z);
				Point3 p(x,y,z);
				ploop.push_back(p);

				// use the closest point's normal as the initial front normal
				MeshGuidanceField::kdtree_type::OrderedIncrementalTraverse ctx(*((MeshGuidanceField*)guidance)->GetKdTree(), p);

				real_type squared_dist;
				int cp = ctx.GetNext(squared_dist);
				Vector3 norm = mesh.verts[cp].normal;
				nloop.push_back(norm);
			}
			fclose(f);

			reverse(ipts);
			ipts.push_back(ploop);
			inorms.push_back(nloop);
		}
	} else {
		// normal boundary extraction
		extract_boundary_loops(mesh, boundaries);
	}

	CreaseExtractor extractor(mesh, (MeshGuidanceField *)guidance,
							  &projector,
							  crease_points, crease_onormals,
							  crease_indices, crease_normals,
							  corner_triangles,
							  sharp,
							  small_crease);
	cerr << crease_indices.size() << endl;
	extractor.Extract();
	cerr << crease_indices.size() << endl;
	//	creases = extractor.Extract();
#endif	

	if (boundaries.size()) {
		cerr << (boundaries.size()) << " boundaries." << endl;
		for (unsigned b=0; b<boundaries.size(); b++) {
			cerr << "Resampling boundary " << b << endl;
			vector<Point3> ploop;
			vector< vector<Vector3> > nloop(1);
			for (unsigned i=0; i<boundaries[b].size(); i++) {
				ploop.push_back(mesh.verts[boundaries[b][i]].point);
				nloop[0].push_back(mesh.verts[boundaries[b][i]].normal);
			}

			vector<Point3> plooprs;
			vector< vector<Vector3> > nlooprs(1);
			guidance->ResampleCurve(ploop, nloop, plooprs, nlooprs);
			
			ipts.push_back(plooprs);
			inorms.push_back(nlooprs[0]);
		}
	} else if (!crease_indices.size() && !ipts.size()) {
		unsigned startvert = (unsigned) 0;
		
		// let the user select the starting location
		if (0 && gui) {
			int win; Point3 p;
			while (!getUISelectPoint(p, win)) { }
			for (unsigned i=1; i<mesh.verts.size(); i++) {
				if (Point3::distance(mesh.verts[i].point, p) < Point3::distance(mesh.verts[startvert].point, p))
					startvert = i;
			}
		}

		cerr<<"starting at point "<<startvert<<endl;


		real_type size = guidance->MaxStepLength(mesh.verts[startvert].point);
		Vector3 udir, vdir;
		PerpVectors(mesh.verts[startvert].normal, udir, vdir);

		ipts.resize(1);
		inorms.resize(1);
		ipts[0].resize(2); inorms[0].resize(2);
		
		ipts[0][0]	 = mesh.verts[startvert].point;
		inorms[0][0] = mesh.verts[startvert].normal;

		projector.ProjectPoint(ipts[0][0] + size * udir, ipts[0][1], inorms[0][1]);
	}


	triangulator = new Triangulator(*controller);
	if (crease_indices.size()) {
	
		triangulator->Go(ipts, inorms, failsafe,
						 crease_points, crease_onormals,
						 crease_indices, crease_normals,
						 corner_triangles);
	} else {
		triangulator->Go(ipts, inorms, failsafe);
	}

	return ret;
}


int do_csg(int argc, char* argv[]) {

	int ret=2;
	int csubdiv=0;
	assert(argc>1 && argv[1][0]!='-');
	char *command = argv[1];
	if (argc>2 && argv[2][0]!='-') {
		csubdiv = atoi(argv[2]);
		ret++;
	}


	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;


	// set the direction of the loops the way we need them
	bool reverse0, reverse1;
	if (stricmp(command, "union")==0) {
		reverse0 = false;
		reverse1 = false;
	} else if (stricmp(command, "int")==0) {
		reverse0 = true;
		reverse1 = true;
	} else if (stricmp(command, "sub")==0) {
		reverse0 = false;
		reverse1 = true;
	} else {
		cerr<<"unknown csg operation: "<<command<<endl;
		return 4;
	}


	vector< vector<Point3> > points0, points1;
	vector< vector<Vector3> > normals0, normals1;
	vector<int> pointsides[2];
	for (int m=0; m<2; m++) {
		pointsides[m].resize(cmeshes[m]->verts.size());
		for (unsigned i=0; i<pointsides[m].size(); i++)
			pointsides[m][i] = 0;	// unknown
	}


	if (cmeshes[0]->num_vertices() && cmeshes[1]->num_vertices()) {
		
		cerr<<"Doing mesh/mesh csg"<<endl;

		cerr<<"getting intersection loops...";
		MeshCSGGuidanceField::GetIntersectionLoops(*cmeshes[0], *cmeshes[1], points0, normals0, normals1, pointsides);
		cerr<<"OK"<<endl;

		if (points0.size() == 0) {
			cerr<<"surfaces don't intersect!"<<endl;
			return 4;
		}

		// figure out which points of the models we need to keep
		cerr<<"flooding point sides...";
		MeshCSGGuidanceField::FloodPointSides(*cmeshes[0], pointsides[0]);
		MeshCSGGuidanceField::FloodPointSides(*cmeshes[1], pointsides[1]);
		cerr<<"OK"<<endl;

		if (!reverse0) {
			for (unsigned i=0; i<pointsides[0].size(); i++)
				pointsides[0][i] = -pointsides[0][i];
		}
		if (reverse1) {
			for (unsigned i=0; i<pointsides[1].size(); i++)
				pointsides[1][i] = -pointsides[1][i];
		}


		guidance = new MeshCSGGuidanceField(csubdiv, *cmeshes[0], *cmeshes[1], pointsides, points0, rho, min_step, max_step, reduction);

		cerr<<"trimming point sides...";
		((MeshCSGGuidanceField*)guidance)->TrimPointSides(*cmeshes[0], pointsides[0]);
		((MeshCSGGuidanceField*)guidance)->TrimPointSides(*cmeshes[1], pointsides[1]);
		cerr<<"OK"<<endl;


		cerr<<"blending curvature...";
		((MeshCSGGuidanceField*)guidance)->BlendGuidance(csg_guidance_blend, cmeshes, pointsides);
		cerr<<"OK"<<endl;


		if (0) {
			dbgClear();
			for (int m=0; m<2; m++) {
				for (unsigned i=0; i<cmeshes[m]->verts.size(); i++) {
					if (pointsides[m][i] == -1) {
						DbgPoints::add(cmeshes[m]->verts[i].point, 0.5, 0, 0);
					} else if (pointsides[m][i] == -2) {
						DbgPoints::add(cmeshes[m]->verts[i].point, 1, 0, 0);
					} else if (pointsides[m][i] == 1) {
						DbgPoints::add(cmeshes[m]->verts[i].point, 0, 0.5, 0);
					} else if (pointsides[m][i] == 2) {
						DbgPoints::add(cmeshes[m]->verts[i].point, 0, 1, 0);
					}
				}
			}
			redrawAndWait(' ');
		}



		// resample the intersection curves
		cerr<<"resampling intersection curves...";
		for (unsigned l=0; l<points0.size(); l++) {
			vector<Point3> op;
			vector< vector<Vector3> > inorms, onorms;
			inorms.push_back(normals0[l]);
			inorms.push_back(normals1[l]);
			guidance->ResampleCurve(points0[l], inorms, op, onorms);

			points0[l] = op;
			normals0[l] = onorms[0];
			normals1[l] = onorms[1];
		}
		cerr<<"OK"<<endl;

		points1 = points0;
		for (unsigned l=0; l<points0.size(); l++) {
			reverse(points1[l]);
			reverse(normals1[l]);
		}

		cerr<<"computing initial fronts...";
		MeshCSGGuidanceField::FrontsFromSides(*cmeshes[0], pointsides[0], points0, normals0);
		MeshCSGGuidanceField::FrontsFromSides(*cmeshes[1], pointsides[1], points1, normals1);
		cerr<<"OK"<<endl;


		// set the direction of the loops the way we need them
		for (unsigned l=0; l<points0.size(); l++) {
			if (reverse0) {
				reverse(points0[l]);
				reverse(normals0[l]);
			}
		}

		for (unsigned l=0; l<points1.size(); l++) {
			if (!reverse1) {
				reverse(points1[l]);
				reverse(normals1[l]);
			}
		}

		controller = new ControllerWrapper(guidance, NULL, GetOutputController());
		triangulator = new Triangulator(*controller);


		MeshProjector projector1(*cmeshes[0], &pointsides[0]);
		controller->SetProjector(&projector1);
		triangulator->SetFlipOutput(reverse0);
		//	hhmout->SetVertexKey("{matid=0}");
		triangulator->InsertSubMesh(*cmeshes[0], pointsides[0], 1);
		//	hhmout->SetVertexKey("{matid=1}");
		triangulator->Go(points0, normals0, failsafe);


		MeshProjector projector2(*cmeshes[1], &pointsides[1]);
		controller->SetProjector(&projector2);
		//		triangulator = new Triangulator(*controller);
		triangulator->SetFlipOutput(reverse1);
		//	hhmout->SetVertexKey("{matid=0}");
		triangulator->InsertSubMesh(*cmeshes[1], pointsides[1], 1);
		//	hhmout->SetVertexKey("{matid=1}");
		triangulator->Go(points1, normals1, failsafe);
 
	} else if (cmeshes[0]->num_vertices() && points.size()) {
		cerr<<"Doing mesh/pointset csg"<<endl;


		MeshProjector mp(*cmeshes[0], NULL);
		SmoothMLSProjector psp(points, adamson);

		compute_pointset_normals(points, psp._projector);

		cerr<<"flip pointset normal orientation? [y/n]"<<endl;
		cerr<<"defaulting n"<<endl;	//for (unsigned i=0; i<points.size(); i++) { points.normal(i).flip(); }
#if 0
		while (1) {
			int key = GetUILastKey();
			if (key=='y' || key=='Y') {
				for (unsigned i=0; i<points.size(); i++)
					points.normal(i).flip();
			}
			break;
		} else if (key=='n' || key=='N') {
			break;
		}
#endif



		cerr<<"computing intersection curves"<<endl;
		vector< vector<Point3> > points0, points1;
		vector< vector<Vector3> > normals0, normals1;
		MeshPSCSGGuidanceField::GetIntersectionLoops(mp, psp, points0, normals0, normals1);
		cerr<<"OK"<<endl;


		guidance = new MeshPSCSGGuidanceField(csubdiv, *cmeshes[0], points, psp._projector, points0, rho, min_step, max_step, reduction, adamson);

		// computing mesh pointsides
		vector<int> pointsides;
		((MeshPSCSGGuidanceField*)guidance)->HackPointSides(*cmeshes[0], pointsides);
		MeshCSGGuidanceField::FixPointSides(*cmeshes[0], pointsides);



		// resample the intersection curves
		cerr<<"resampling intersection curves...";
		for (unsigned l=0; l<points0.size(); l++) {
			vector<Point3> op;
			vector< vector<Vector3> > inorms, onorms;
			inorms.push_back(normals0[l]);
			inorms.push_back(normals1[l]);
			guidance->ResampleCurve(points0[l], inorms, op, onorms);

			points0[l] = op;
			normals0[l] = onorms[0];
			normals1[l] = onorms[1];
		}
		cerr<<"OK"<<endl;

		points1 = points0;
		for (unsigned l=0; l<points0.size(); l++) {
			reverse(points0[l]);
			reverse(normals0[l]);
		}

		cerr<<"computing initial fronts...";
		MeshCSGGuidanceField::FrontsFromSides(*cmeshes[0], pointsides, points0, normals0);
		cerr<<"OK"<<endl;



		// set the direction of the loops the way we need them
		for (unsigned l=0; l<points0.size(); l++) {
			if (reverse0) {
				reverse(points0[l]);
				reverse(normals0[l]);
			}
		}

		for (unsigned l=0; l<points1.size(); l++) {
			if (reverse1) {
				reverse(points1[l]);
				reverse(normals1[l]);
			}
		}


		controller = new ControllerWrapper(guidance, NULL, GetOutputController());
		triangulator = new Triangulator(*controller);


		controller->SetProjector(&mp);
		triangulator->SetFlipOutput(reverse0);
		//	hhmout->SetVertexKey("{matid=0}");
		triangulator->InsertSubMesh(*cmeshes[0], pointsides, 1);
		//	hhmout->SetVertexKey("{matid=1}");
		triangulator->Go(points0, normals0, failsafe);

		controller->SetProjector(&psp);
		triangulator->SetFlipOutput(reverse1);
		//	hhmout->SetVertexKey("{matid=2}");
		//		triangulator->InsertSubMesh(*cmeshes[1], pointsides[1], 1);
		triangulator->Go(points1, normals1, failsafe);

	}

	return ret;
}


int do_tri_vol(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	real_type isoval = atof(argv[1]) * iso_scale;
	bool bspline=false;
	if (argc>2 && !stricmp(argv[2], "bspline"))
		bspline=true;

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;



	IsoSurfaceProjector projector(*cvolume, isoval);

	critical_section->enter();
	volume->SetSplineCoeffs(bspline);
	critical_section->leave();


	if (is_segmentation) {
		guidance = new IsoSurfaceGuidanceField(projector, *(MultiMaterialVolume*)volume, rho, min_step, max_step, reduction);
	} else {
		guidance = new IsoSurfaceGuidanceField(projector, *volume, rho, min_step, max_step, reduction);
	}

	controller = new ControllerWrapper(guidance, &projector, GetOutputController());



	triangulator = new Triangulator(*controller);
	triangulator->SetFlipOutput(true);

	vector< vector<Point3> > ipts;
	vector< vector<Vector3> > inorms;
	triangulator->Go(ipts, inorms, failsafe);

	cerr<<volume->GetNumBlockLoads()<<" block loads"<<endl;

	if (bspline)
		return 3;
	else
		return 2;
}


int do_make_gf(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-' && argv[2][0]!='-');
	int block = atof(argv[1]);
	real_type isoval = atof(argv[2]) * iso_scale;

	bool bspline=false;
	if (argc>3 && !stricmp(argv[3], "bspline"))
		bspline=true;

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;


	IsoSurfaceProjector projector(*cvolume, isoval);

	critical_section->enter();
	volume->SetSplineCoeffs(bspline);
	critical_section->leave();


	// reuse_guidance = false; // always make new this way
	// reuse_seeds = false; // always make new this way
	guidance = new IsoSurfaceGuidanceField(projector, *volume, rho, min_step, max_step, reduction, block);


	if (bspline)
		return 4;
	else
		return 3;
}



int do_tri_tet(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	real_type isoval = atof(argv[1]) * iso_scale;

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;



	critical_section->enter();
	vector<real_type> scalars;
	ctetmesh->GetShell(tetshell, scalars);
	for (int i=0; i<0; i++) {
		tetshell.LoopSubdivide(true);
	}
	tetshell.BuildStrips();
	critical_section->leave();


	redrawAndWait(' ');

	TetMeshProjector *projector = (nielson!=0 ?
								   (TetMeshProjector*) new TetMeshProjectorNielson(*ctetmesh, isoval) :
								   (TetMeshProjector*) new TetMeshProjectorMLS(*ctetmesh, isoval));


	guidance = new TetMeshGuidanceField(*projector, rho, min_step, max_step, reduction);
	controller = new ControllerWrapper(guidance, projector, GetOutputController());


	allow_outside = true;


	vector< vector<Point3> > ipts;
	vector< vector<Vector3> > inorms;

	if (nielson == 0) {
		cerr<<"evaluating at shell points"<<endl;
		scalars.resize(ctetshell->verts.size());
		for (unsigned v=0; v<ctetshell->verts.size(); v++) {
			projector->EvalAtPoint(ctetshell->verts[v].point, scalars[v], NULL);
			//		cerr<<scalars[v]<<endl;
		}
	}

	redrawAndWait(' ');

	cerr<<"doing marching triangles"<<endl;
	MeshProjector shell_projector(*ctetshell);
	vector< vector<Point3> > ploops;
	MarchingTris(*ctetshell, *projector, scalars, isoval, ploops);


	if (ploops.size()) {

		vector< vector<Vector3> > nloops(ploops.size());
		for (unsigned l=0; l<ploops.size(); l++) {
			nloops[l].resize(ploops[l].size());


			if (nielson != 0) {
				// move the points so we're not trying to evaluate on teh edges of the tet!
				vector<Point3> pl(ploops[l].size());
				for (unsigned p=0; p<ploops[l].size(); p++) {
					pl[p] = Point3::inbetween(ploops[l][p], ploops[l][(p+1)%ploops[l].size()], 0.01);
				}
				ploops[l] = pl;
			}


			// project the points onto the isosurface / get the normals
			for (unsigned p=0; p<ploops[l].size(); p++) {
				DbgPoints::add(ploops[l][p], 1,0,0);
				/*
				  if (projector.ProjectPoint(ploops[l][p], ploops[l][p], nloops[l][p]) != PROJECT_SUCCESS) {
				  cerr<<"couldn't project boundary point!"<<endl; return 2;
				  }
				*/

				projector->NormalAtPoint(ploops[l][p], nloops[l][p]);
				DbgPoints::add(ploops[l][p], 0,1,0);

				vector<Point3> line(2);
				line[0] = ploops[l][p];
				line[1] = line[0] + (real_type)0.1 * nloops[l][p];
				DbgPLines::add(line, 0);

			}
			redrawAndWait(' ');


			// resample
			while (1) {
				vector< vector<Vector3> > nloop2(1);	nloop2[0] = nloops[l];
				vector<Point3> plooprs;
				vector< vector<Vector3> > nlooprs(1);
				guidance->ResampleCurve(ploops[l], nloop2, plooprs, nlooprs);


				cerr<<"resampled loop "<<l<<endl;
				dbgClear();
				for (unsigned i=0; i<plooprs.size(); i++) {
					DbgPoints::add(plooprs[i], 1,0,0);
				}

				DbgPLines::add(plooprs, 1);
				redrawAndWait(' ');


				vector<Point3> oplooprs = plooprs;
				for (unsigned p=0; p<plooprs.size(); p++) {

					if (p>0) {
						Vector3 lastmove = plooprs[p-1] - oplooprs[p-1];
						//				plooprs[p] += lastmove;
					}

					Vector3 meshnorm;
					Vector3 oldnorm = nlooprs[0][p];
					if (!ProjectToIntersection(shell_projector, *projector, plooprs[p], meshnorm, nlooprs[0][p], 1e-4)) {
						nlooprs[0][p] = oldnorm;
					}
				}

				if (plooprs.size() <= ploops[l].size()) {

					for (int i=0; i<plooprs.size(); i++) {
						real_type d1 = Point3::distance(plooprs[i], plooprs[(i+plooprs.size()-1)%plooprs.size()]);
						real_type d2 = Point3::distance(plooprs[i], plooprs[(i+1)%plooprs.size()]);
						real_type d3 = Point3::distance(plooprs[(i+plooprs.size()-1)%plooprs.size()], plooprs[(i+1)%plooprs.size()]);

						if (1.1*d3<d1 || 1.1*d3<d2 ||
							d1 < d2*0.2) {
							plooprs.erase(plooprs.begin()+i);
							nlooprs[0].erase(nlooprs[0].begin()+i);
							i--;
						}
					
					}

					ploops[l] = plooprs;
					nloops[l] = nlooprs[0];
					break;
				}


				ploops[l] = plooprs;
				nloops[l] = nlooprs[0];
			}


			if (ploops[l].size() > 2) {
				cerr<<"adding loop "<<l<<endl;

				dbgClear();
				for (unsigned i=0; i<ploops[l].size(); i++) {
					DbgPoints::add(ploops[l][i], 1,0,0);
				}
				redrawAndWait(' ');

				ipts.push_back(ploops[l]);
				inorms.push_back(nloops[l]);
			}
		}

		// check for some weirdness that happens on the spx
		vector<int> gaps;
		for (unsigned l=0; l<ipts.size(); l++) {
			int n = ipts[l].size();
			vector<Point3> &loop = ipts[l];


			real_type longest=0;
			int li = 0;
			for (int i=0; i<n; i++) {
				real_type d = Point3::distance(loop[(i+1)%n], loop[i]);
				if (d > longest) {
					longest=d;
					li = i;
				}
			}

			real_type longest2=0;
			int li2 = 0;
			for (int i=0; i<n; i++) {
				if (i==li) continue;
				real_type d = Point3::distance(loop[(i+1)%n], loop[i]);
				if (d > longest2) {
					longest2=d;
					li2 = i;
				}
			}


			real_type ratio = Point3::distance(loop[(li+1)%n], loop[li]) /
				std::max(Point3::distance(loop[(li+n-1)%n], loop[li]),
						 Point3::distance(loop[(li+1)%n], loop[(li+2)%n]));

			if (longest > 1.5*longest2 &&
				ratio > 1.5) {
				gaps.push_back(l);
				gaps.push_back(li);
			}
		}

		if (gaps.size() == 4) {

			cerr<<"merging loops!"<<endl;

			int l1 = gaps[0];
			int li1 = gaps[1];
			int l2 = gaps[2];
			int li2 = gaps[3];

			vector<Point3> prot;
			vector<Vector3> nrot;

			prot.insert(prot.end(), ipts[l1].begin()+li1+1, ipts[l1].end());
			prot.insert(prot.end(), ipts[l1].begin(), ipts[l1].begin()+li1+1);

			nrot.insert(nrot.end(), inorms[l1].begin()+li1+1, inorms[l1].end());
			nrot.insert(nrot.end(), inorms[l1].begin(), inorms[l1].begin()+li1+1);

			ipts[l1] = prot;
			inorms[l1] = nrot;


			prot.clear();
			nrot.clear();
			prot.insert(prot.end(), ipts[l2].begin()+li2+1, ipts[l2].end());
			prot.insert(prot.end(), ipts[l2].begin(), ipts[l2].begin()+li2+1);

			nrot.insert(nrot.end(), inorms[l2].begin()+li2+1, inorms[l2].end());
			nrot.insert(nrot.end(), inorms[l2].begin(), inorms[l2].begin()+li2+1);

			ipts[l2] = prot;
			inorms[l2] = nrot;

			
			ipts[l1].insert(ipts[l1].end(), ipts[l2].begin()+1, ipts[l2].end()-1);
			inorms[l1].insert(inorms[l1].end(), inorms[l2].begin()+1, inorms[l2].end()-1);

			swap(ipts[l2], ipts.back());
			swap(inorms[l2], inorms.back());

			ipts.pop_back();
			inorms.pop_back();


			vector< vector<Vector3> > nloop2(1);	nloop2[0] = inorms[l1];
			vector<Point3> plooprs;
			vector< vector<Vector3> > nlooprs(1);
			guidance->ResampleCurve(ipts[l1], nloop2, plooprs, nlooprs);
			ipts[l1] = plooprs;
			inorms[l1]=  nlooprs[0];

		}

	} else {


		int startvert = 0;

		Point3 startpoint = ((TetMeshGuidanceField*)guidance)->PointLocation(startvert);
		Vector3 startnorm;
		if (!projector->NormalAtPoint(startpoint, startnorm)) {
			cerr<<"couldn't get normal"<<endl;
			exit(0);
		}
		real_type size = guidance->MaxStepLength(startpoint);

		Vector3 udir, vdir;
		PerpVectors(startnorm, udir, vdir);

		ipts.resize(1);
		inorms.resize(1);
		ipts[0].resize(2); inorms[0].resize(2);

		ipts[0][0]	 = startpoint;
		inorms[0][0] = startnorm;
		projector->ProjectPoint(ipts[0][0] + size * udir, ipts[0][1], inorms[0][1]);;
	}


	allow_outside = false;


	triangulator = new Triangulator(*controller);
	triangulator->SetFlipOutput(true);
	triangulator->Go(ipts, inorms, failsafe);

	delete projector;

	return 2;
}


int do_tet_smooth(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	int passes = atoi(argv[1]);

	critical_section->enter();
	tetmesh.SmoothScalars(passes);
	critical_section->leave();

	return 2;
}


void do_volsize() {
	int n=0;
	for (int x=0; x<volume->GetDim(0); x++) {
		for (int y=0; y<volume->GetDim(1); y++) {
			for (int z=0; z<volume->GetDim(2); z++) {
				float v = volume->GetValue(x,y,z);
				if (v != 0)
					n++;
			}
		}
	}
	cout<<n<<endl;
}


int do_segmentation(int argc, char* argv[]) {

	if (argc>1 && argv[1][0] != '-') {

		MultiMaterialVolume *rv = new MultiMaterialVolume(*volume);
		delete volume;
		volume = rv->matProbs[atoi(argv[1])];
		cvolume= volume;

		return 2;
	} else {

		MultiMaterialVolume *rv = new MultiMaterialVolume(*volume);

		delete volume;

		volume = rv;
		cvolume= rv;

		is_segmentation = true;

		return 1;
	}
}

int do_vol_smooth(int argc, char* argv[]) {

	assert(argc>1 && argv[1][0]!='-');
	int width = atoi(argv[1]);

	RegularVolume<float> *rv = new RegularVolume<float>();
	rv->Smoothed(*volume, width);

	delete volume;

	volume = rv;
	cvolume= rv;


	return 2;
}


int do_marchingcubes(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	bool bspline=false;
	if (argc>2 && !stricmp(argv[2], "bspline"))
		bspline=true;

	critical_section->enter();
	volume->SetSplineCoeffs(bspline);
	critical_section->leave();

	real_type isoval = 0;
	if (!stricmp(argv[1], "random")) {
		srand(time(0));
		while (isoval == 0) {
			isoval = volume->EvalAtPoint(Point3((rand()/(float)RAND_MAX)*volume->GetAspect(0)*volume->GetDim(0),
												(rand()/(float)RAND_MAX)*volume->GetAspect(1)*volume->GetDim(1),
												(rand()/(float)RAND_MAX)*volume->GetAspect(2)*volume->GetDim(2)) +
										 volume->GetTranslation());
		}
		//		cerr<<isoval<<endl;
	} else {
		isoval = iso_scale*atof(argv[1]);
	}


  
	extern real_type disoval;
	disoval=isoval;
	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;

	OutputController *oc = NULL;

	if (mc_stats == 0) {
		oc = GetOutputController(false);
		controller = new ControllerWrapper(NULL, NULL, oc);	// for the output
	} else if (mc_stats == 2) {
		if (gui) {
			oc = gui;
			controller = new ControllerWrapper(NULL, NULL, oc);	// for the output
		} else {
			oc = NULL;
			controller = NULL;
		}
	}


	if (mc_stats == 0) {
		MarchingCubes(*cvolume, *controller, isoval, true);
	} else {
		real_type cells, wcells, tris, wtris, area, warea;
		MarchingCubesStats(*cvolume, controller, isoval,
						   cells, wcells, tris, wtris, area, warea);

		fileappend(".uarea.txt", area);
		fileappend(".warea.txt", warea);
		fileappend(".utris.txt", tris);
		fileappend(".wtris.txt", wtris);
		fileappend(".ucells.txt", cells);
		fileappend(".wcells.txt", wcells);

	}

	if (oc)
		oc->Finish();


	if (bspline)
		return 3;
	else
		return 2;
}


int do_marchingtets(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	real_type isoval = atof(argv[1]) * iso_scale;


	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;

	// never use edgeflipper on mt
	controller = new ControllerWrapper(NULL, NULL, GetOutputController(false, true));	// for the output
	MarchingTets(*ctetmesh, *controller, isoval);

	return 2;
}


int do_tri_mls(int argc, char* argv[]) {
	//	  assert(argc>1 && argv[1][0]!='-');

	int ret=1;

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;

	//	  read_obj(argv[1], points);


	SmoothMLSProjector projector(points, adamson);
	projector.SetRadiusFactor(radius_factor);

	// compute the normals for the adamson projection
	if (adamson==1 && !points.has_normals()) {
		compute_pointset_normals(points, projector._projector);
	}

	guidance = new SmoothMLSGuidanceField(projector._projector, rho, min_step, max_step, reduction, adamson, gffile);
	controller = new ControllerWrapper(guidance, &projector, GetOutputController());


	vector< vector<Point3> > ipts;
	vector< vector<Vector3> > inorms;


	if (argc >= 2 && endswith(argv[1], ".flf")) {
		ret++;

		// read Joel/Linh's feature file and use it as initial fronts
		vector<InfoPoint> pts;
		vector< vector<InfoIndex> > loops;
		FLF_read(argv[1], pts, loops);


		for (unsigned i=0; i<loops.size(); i++) {

			int start_corner = 0;
			for (unsigned j=0; j<loops[i].size(); j++) {
				if (pts[loops[i][j].index_].isCorner_) {
					start_corner = j; break;
				}
			}

			ipts.push_back(vector<Point3>());
			inorms.push_back(vector<Vector3>());


			int first_corner = start_corner;

			while (1) {

				vector<Point3> points;
				vector< vector<Vector3> > normals(1);
				vector<Point3> rpoints;
				vector< vector<Vector3> > rnormals(1);


				int end_corner = start_corner;
				while (1) {
					points.push_back(Point3(pts[loops[i][end_corner].index_].x_,
											pts[loops[i][end_corner].index_].y_,
											pts[loops[i][end_corner].index_].z_));
					normals[0].push_back(Vector3(loops[i][end_corner].x_,
												 loops[i][end_corner].y_,
												 loops[i][end_corner].z_));

					end_corner = (end_corner+1) % loops[i].size();
			
					if (end_corner == start_corner) break;
					if (pts[loops[i][end_corner].index_].isCorner_) break;
				}
				points.push_back(Point3(pts[loops[i][end_corner].index_].x_,
										pts[loops[i][end_corner].index_].y_,
										pts[loops[i][end_corner].index_].z_));
				normals[0].push_back(Vector3(loops[i][end_corner].x_,
											 loops[i][end_corner].y_,
											 loops[i][end_corner].z_));


				// see if we should reverse the ordering just for resampling
				bool doreverse = (loops[i][first_corner].index_ > loops[i][end_corner].index_);
				if (doreverse) {
					reverse(points);
					reverse(normals[0]);
				}

				guidance->ResampleCurve(points, normals, rpoints, rnormals, true);

				// check if we have to flip the normals back over
				if (rnormals[0][0].dot(normals[0][0]) < 0) {
					cerr<<"flipping normals of loop "<<i<<endl;
					for (unsigned j=0; j<rnormals[0].size(); j++) {
						rnormals[0][j] *= -1;
					}
				}

				// reverse the order back
				if (doreverse) {
					reverse(rpoints);
					reverse(rnormals[0]);
				}

				for (unsigned j=0; j<rpoints.size()-1; j++) {
					ipts.back().push_back(rpoints[j]);
					inorms.back().push_back(rnormals[0][j]);
				}
		
				start_corner = end_corner;
				if (start_corner == first_corner) break;
			}


			// now make sure they're oriented right
			int lproj=0;
			int rproj=0;
			real_type ldist=0;
			real_type rdist=0;
			for (int i=0; i<ipts.back().size(); i++) {

				Vector3 dir = ipts.back()[(i+1)%ipts.back().size()] - ipts.back()[i];
				Vector3 norm = inorms.back()[i];
				Vector3 left = norm.cross(dir);

				Point3 fp, tp;
				Vector3 tn;

				fp = ipts.back()[i] + left;
				DbgPoints::add(fp, 0,1,0);
				if (projector.ProjectPoint(fp, tp, tn) == PROJECT_SUCCESS) {
					lproj++;
					ldist = Point3::distance(fp, tp);
				}

				fp = ipts.back()[i] - left;
				DbgPoints::add(fp, 0,0,1);
				if (projector.ProjectPoint(fp, tp, tn) == PROJECT_SUCCESS) {
					rproj++;
					rdist = Point3::distance(fp, tp);
				}
			}

			cerr<<rproj<<" "<<(rdist/rproj)<<endl;
			cerr<<lproj<<" "<<(ldist/lproj)<<endl;

			if (rdist/rproj < ldist/lproj) {
				cerr<<"reversing"<<endl;
				reverse(ipts.back());
				reverse(inorms.back());
			} else {
				cerr<<"not reversing"<<endl;
			}
		}

		redrawAndWait(' ');

	} else {
		/*! Compute an initial front
		 * This is simply a random point + another adjacent point
		 * which forms an edge
		 */
		//	Point3 x0 = points.vertex((int)(gtb::nrran1f()*points.size()));
		Point3 x0 = points.vertex(0);

		Vector3 n0;
		projector.ProjectPoint(x0, x0, n0);
		real_type len = guidance->MaxStepLength(x0);
		Vector3 udir, vdir;
		PerpVectors(n0, udir, vdir);
		Point3 x1 = x0 + udir*len;
		Vector3 n1;
		projector.ProjectPoint(x1, x1, n1);
		if (n1.dot(n0) < 0.0f) n1.flip();

		ipts.resize(1);
		inorms.resize(1);
		ipts[0].resize(2); inorms[0].resize(2);
		ipts[0][0] = x0;
		ipts[0][1] = x1;
		inorms[0][0] = n0;
		inorms[0][1] = n1;
	}

	triangulator = new Triangulator(*controller);
	triangulator->Go(ipts, inorms, failsafe);

	return ret;
}



int do_tri(int argc, char* argv[]) {

	if (cmeshes[0]->verts.size())
		return do_tri_mesh(argc, argv);
	else if (cvolume)
		return do_tri_vol(argc, argv);
	else if (points.size())
		return do_tri_mls(argc, argv);
	else if (tetmesh.verts.size())
		return do_tri_tet(argc, argv);

	cerr<<"nothing loaded to triangulate!"<<endl;
	return 1;
}


int do_subdiv(int argc, char* argv[]) {

	// since we need to modify the mesh
	critical_section->enter();

	assert(argc>1 && argv[1][0]!='-');
	int num = atoi(argv[1]);

	for (int i=0; i<num; i++) {
		meshes[0].LoopSubdivide();
		meshes[1].LoopSubdivide();
	}

	if (meshes[0].verts.size())
		meshes[0].BuildStrips(100);
	if (meshes[1].verts.size())
		meshes[1].BuildStrips(100);


	critical_section->leave();

	return 2;
}


void do_flip() {

	critical_section->enter();

	for (int f=0; f<(int)meshes[0].faces.size(); f++) {
		std::swap(meshes[0].faces[f].verts[1], meshes[0].faces[f].verts[2]);
		std::swap(meshes[0].faces[f].nbrs[1], meshes[0].faces[f].nbrs[2]);
		std::swap(meshes[0].faces[f].nbrs[2], meshes[0].faces[f].nbrs[0]);
		std::swap(meshes[0].faces[f].nbrs[0], meshes[0].faces[f].nbrs[1]);
	}

	for (int v=0; v<(int)meshes[0].verts.size(); v++) {
		meshes[0].verts[v].normal.flip();
	}

	meshes[0].BuildStrips(100);

	critical_section->leave();
}


void do_flip_pts() {

	if (gui) {
		SmoothMLSProjector psp(points, adamson);
		psp.SetRadiusFactor(radius_factor);

		cerr<<"ok"<<endl;

		while (1) {
			Point3 p;
			int w;

			bool done = false;
			while (1) {
				int key = GetUILastKey();
				if (key == ' ') {
					if (getUISelectPoint(p, w)) {
						break;
					}
				} else if (key == 'd' || key == 'D') {
					cerr<<"finished"<<endl;
					done=true; break;
				} else if (key == 'v' || key == 'V') {
					cerr<<"flipping by ave"<<endl;
					vector<Vector3> norms(points.size());
					for (int i=0; i<norms.size(); i++) {
						surfelset_view nbhd(points);
						psp._projector.extract(points.vertex(i), ((real_type)1)*psp._projector.point_radius(points.vertex(i)), nbhd);

						Vector3 avenorm(0,0,0);
						for (int j=0; j<nbhd.size(); j++) {
							avenorm += points.normal(nbhd.get_index(j));
						}
						avenorm.normalize();

						if (avenorm.dot(points.normal(i)) > 0.3)
							norms[i] = points.normal(i);
						else if (avenorm.dot(points.normal(i)) < -0.3)
							norms[i] = -points.normal(i);
						else
							norms[i] = avenorm;
						
					}

					for (int i=0; i<norms.size(); i++) {
						points.normal(i) = norms[i];
					}

					cerr<<"done"<<endl;
				}
			}
			if (done)
				break;

			int closest_i = 0;
			real_type closest_d = INFINITY;
			for (int i=0; i<points.size(); i++) {
				real_type d = Point3::distance(points.vertex(i), p);
				if (d < closest_d) {
					closest_d = d;
					closest_i = i;
				}
			}

			cerr<<"starting flip"<<endl;
			points.normal(closest_i) = -points.normal(closest_i);
			orient_normals(closest_i, points, psp._projector, false);
			cerr<<"done"<<endl;
		}
		
	} else {
		for (int i=0; i<points.normals().size(); i++) {
			points.normal(i) = -points.normal(i);
		}
	}
}



// interactively let the user flip the orientations of the connected components
void do_flip_components() { 

	if (!gui) {
		cerr<<"gui needed for interactive component flipping!"<<endl;
		return;
	}

	vector<int> done(meshes[0].faces.size(), 0);

	while (1) {
		vector<int> toprocess;

		for (unsigned i=0; i<meshes[0].faces.size(); i++) {
			if (!done[i]) {
				toprocess.push_back(i);
				break;
			}
		}

		if (toprocess.empty())
			break;

		vector<int> component;

		while (!toprocess.empty()) {

			int next = toprocess.back();
			toprocess.resize(toprocess.size()-1);

			if (done[next])
				continue;

			done[next] = true;
			component.push_back(next);

			for (int i=0; i<3; i++) {
				if (meshes[0].faces[next].nbrs[i] >= 0 &&
					!done[meshes[0].faces[next].nbrs[i]])
					toprocess.push_back(meshes[0].faces[next].nbrs[i]);
			}

		}


		dbgClear();
		for (unsigned i=0; i<component.size(); i++) {
			DbgPoly::add(meshes[0].verts[meshes[0].faces[component[i]].verts[0]].point, Point3(1,0,0),
						 meshes[0].verts[meshes[0].faces[component[i]].verts[1]].point, Point3(1,0,0),
						 meshes[0].verts[meshes[0].faces[component[i]].verts[2]].point, Point3(1,0,0));

		}
		gui->Redraw();


		bool stop = false;
		while (!stop) {
			gui->WaitForKey();
			int key = gui->LastKey();

			switch (key) {
			case 'f':
			case 'F':
				{

					std::set<int> verts;
					for (unsigned i=0; i<component.size(); i++) {
						for (int j=0; j<3; j++) {
							verts.insert(meshes[0].faces[component[i]].verts[j]);
						}
					}

					critical_section->enter();
					for (unsigned i=0; i<component.size(); i++) {
						int f = component[i];
						std::swap(meshes[0].faces[f].verts[1], meshes[0].faces[f].verts[2]);
						std::swap(meshes[0].faces[f].nbrs[1], meshes[0].faces[f].nbrs[2]);
						std::swap(meshes[0].faces[f].nbrs[2], meshes[0].faces[f].nbrs[0]);
						std::swap(meshes[0].faces[f].nbrs[0], meshes[0].faces[f].nbrs[1]);
					}
			
			
					for (std::set<int>::iterator v=verts.begin(); v!=verts.end(); ++v) {
						meshes[0].verts[*v].normal.flip();
					}
			
					meshes[0].BuildStrips(100);
					critical_section->leave();
					gui->Redraw();
				}
				break;

			case ' ':
				stop=true;
				break;
			}

		}
	}
}



void do_trans() {

	if (!gui) {
		cerr<<"gui must be used to transform meshes!"<<endl;
		exit(-1);
	}

	gui->MeshEditMode(meshes);
}


int do_save_pts(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	// gcc 4.1 can't figure out where this function is for some reason
	//write_obj(argv[1], points);
	//	cerr<<"CAN'T DO SAVE_PTS!"<<endl;

	FILE *f = fopen(argv[1], "w");
	for (int i=0; i<points.size(); i++) {
		Point3 p = points.vertex(i);
		fprintf(f, "v %g %g %g\n", p[0], p[1], p[2]);

		if (points.has_normals()) {
			Vector3 n = points.normal(i);
			fprintf(f, "vn %g %g %g\n", n[0], n[1], n[2]);
		}
	}

	fclose(f);

	return 2;
}


void do_save_r3d() {
	for (int v=0; v<cmeshes[0]->verts.size(); v++) {

		if (cmeshes[0]->verts[v].normal.length() < 0.5) {
			continue;
		}

		cout<<"v "<<cmeshes[0]->verts[v].point[0]<<" "
			<<cmeshes[0]->verts[v].point[1]<<" "
			<<cmeshes[0]->verts[v].point[2]<<" "
			<<"\nvn "<<cmeshes[0]->verts[v].normal[0]<<" "
			<<cmeshes[0]->verts[v].normal[1]<<" "
			<<cmeshes[0]->verts[v].normal[2]<<endl;
	}
}


int do_save_mesh0(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	cmeshes[0]->Write(argv[1]);
	return 2;
}

int do_save_mesh1(int argc, char* argv[]) {
	assert(argc>1 && argv[1][0]!='-');
	cmeshes[1]->Write(argv[1]);
	return 2;
}

template <typename T>
bool contains(const vector<T> &v, const T &x) {
	return (find(v.begin(), v.end(), x) != v.end());
}

int do_multi_material(int argc, char* argv[]) {

	int used_args = 1;
	bool bspline=false;
	if (argc>1 && !stricmp(argv[1], "bspline")) {
		bspline=true;
		used_args++;
	}

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;


	vector<TriangleMesh*> meshes;
	vector<RegularVolume<float>* > volumes;
	for (int i=used_args; i<argc; i+=2) {
		if (argv[i][0] == '-') break;
		meshes.push_back(new TriangleMesh());
		meshes.back()->Read(argv[i]);
		volumes.push_back((RegularVolume<float>*)VolumeRead(argv[i+1]));
		used_args+=2;
	}

	vector< vector<Point3> > loops;
	vector< vector<int> > loop_mats;

	for (int mat0=0; mat0<meshes.size()-1; mat0++) {
		for (int mat1=mat0+1; mat1<meshes.size(); mat1++) {

			cerr<<mat0<<" "<<mat1<<endl;
			
			cerr<<"kd tree"<<endl;
			vector<Point3> m2verts(meshes[mat1]->verts.size());
			for (int i=0; i<m2verts.size(); i++)
				m2verts[i] = meshes[mat1]->verts[i].point;
			
			typedef gtb::KDTree<int, real_type, GetPointVector<Point3> > kdtree_type;
			kdtree_type kd(10, Box3::bounding_box(m2verts), GetPointVector<Point3>(m2verts));
			for (int i=0; i<m2verts.size(); i++)
				kd.Insert(i);
			kd.MakeTree();

			// mark all vertices that are shared between these two meshes
			cerr<<"marking interface"<<endl;
			bool have_interface = false;
			vector<bool> interface(meshes[mat0]->verts.size(), false);
			for (int v0=0; v0<meshes[mat0]->verts.size(); v0++) {
				kdtree_type::OrderedIncrementalTraverse ot(kd, meshes[mat0]->verts[v0].point);
				real_type squared_dist;
				int v1 = ot.GetNext(squared_dist);
				{
				//				for (int v1=0; v1<meshes[mat1]->verts.size(); v1++) {
					if (meshes[mat0]->verts[v0].someface >= 0 &&
						meshes[mat0]->verts[v0].point == meshes[mat1]->verts[v1].point) {
						interface[v0] = true;
						have_interface = true;
					}
				}
			}

			// pull out the boundary regions of the marked area
			cerr<<"finding boundaries"<<endl;
			bool have_boundary = false;
			vector<bool> done(meshes[mat0]->verts.size(), false);
			for (int v=0; v<meshes[mat0]->verts.size(); v++) {
				if (done[v]) continue;
				if (!interface[v]) {
					done[v] = true;
					continue;
				}

				// see if we're on the boundary
				bool have_marked = false;
				bool have_unmarked = false;
				
				for (TriangleMesh::VertexVertexIteratorI vi(*meshes[mat0], v); !vi.done(); ++vi) {
					int nv = *vi;
					if (interface[nv])
						have_marked = true;
					else
						have_unmarked = true;
				}

				if (!have_marked || !have_unmarked) {
					done[v] = true;
					continue;
				}


				// start a new loop from here
				have_boundary = true;
				cerr<<"making loop"<<endl;
				vector<Point3> loop;
				loop.push_back(meshes[mat0]->verts[v].point);
				done[v] = true;

				int lv = v;
				
				while(1) {

					vector<int> ring;
					for (TriangleMesh::VertexVertexIteratorI vi(*meshes[mat0], lv); !vi.done(); ++vi) {
						ring.push_back(*vi);
					}

					int next = -1;
					for (int i=0; i<ring.size(); i++) {
						if (!interface[ring[i]] && interface[ring[(i+1)%ring.size()]]) {
							next = ring[(i+1)%ring.size()];
						}
					}

					if (next == -1) {
						cerr<<"wtf"<<endl;
					}

					if (v == next)
						break;

					lv = next;

					loop.push_back(meshes[mat0]->verts[next].point);
					done[next] = true;

					DbgPoints::add(meshes[mat0]->verts[next].point, 0,1,0);
				}


				bool loop_found = false;
				for (int l=0; l<loops.size(); l++) {
					if (contains(loops[l], loop.front())) {
						loop_found = true;
						if (!contains(loop_mats[l], mat0))
							loop_mats[l].push_back(mat0);
						if (!contains(loop_mats[l], mat1))
							loop_mats[l].push_back(mat1);
					}
				}

				if (!loop_found) {
					loops.push_back(loop);
					loop_mats.push_back(vector<int>());
					loop_mats.back().push_back(mat0);
					loop_mats.back().push_back(mat1);
				}
			}


			// go back through and flood anything that has an interface boundary
			for (int v=0; v<meshes[mat0]->verts.size(); v++) {
				if (!interface[v])
					continue;

				bool have_marked = false;
				bool have_unmarked = false;
				
				for (TriangleMesh::VertexVertexIteratorI vi(*meshes[mat0], v); !vi.done(); ++vi) {
					int nv = *vi;
					if (interface[nv])
						have_marked = true;
					else
						have_unmarked = true;
				}

				if (!have_marked || !have_unmarked) {
					continue;
				}


				vector<int> toflood;
				toflood.push_back(v);
				interface[v] = false;

				while (!toflood.empty()) {
					int fv = toflood.back();
					toflood.pop_back();
					
					for (TriangleMesh::VertexVertexIteratorI vi(*meshes[mat0], fv); !vi.done(); ++vi) {
						int nv = *vi;
						if (interface[nv]) {
							toflood.push_back(nv);
							interface[nv] = false;
						}
					}
				}
			}


			// go back through again - everything else is connected components
			for (int v=0; v<meshes[mat0]->verts.size(); v++) {
				if (!interface[v])
					continue;

				int count = 1;

				vector<int> toflood;
				toflood.push_back(v);
				interface[v] = false;

				while (!toflood.empty()) {
					int fv = toflood.back();
					toflood.pop_back();
					
					for (TriangleMesh::VertexVertexIteratorI vi(*meshes[mat0], fv); !vi.done(); ++vi) {
						int nv = *vi;
						if (interface[nv]) {
							count++;
							toflood.push_back(nv);
							interface[nv] = false;
						}
					}
				}

				cerr<<count<<endl;
				if (count > 1) {
					loops.push_back(vector<Point3>());
					loops.back().push_back(meshes[mat0]->verts[v].point);

					loop_mats.push_back(vector<int>());
					loop_mats.back().push_back(mat0);
					loop_mats.back().push_back(mat1);
					DbgPoints::add(loops.back().back(), 1,0,0);
				}
			}
		}
	}


	for (int l=0; l<loops.size(); l++) {
		DbgPLines::add(loops[l], 0);
		cerr<<"loop "<<l<<endl;
		for (int l2=0; l2<loop_mats[l].size(); l2++) {
			cerr<<loop_mats[l][l2]<<" ";
		}
		cerr<<endl;
	}
	redrawAndWait(' ');


	// run the triangulation for each material interface
	// we don't do this from tri_vol since we have extra info we need (loops)

	MultiMaterialVolume *mmvol = new MultiMaterialVolume(volumes);
	mmvol->SetSplineCoeffs(bspline);

	IsoSurfaceProjector projector(*mmvol, 0);
	guidance = new IsoSurfaceGuidanceField(projector, *mmvol, loops, loop_mats, rho, min_step, max_step, reduction);
	controller = new ControllerWrapper(guidance, &projector, GetOutputController());

	int num_faces = 0;
	int num_verts = 0;

	for (int mat0=0; mat0<meshes.size()-1; mat0++) {
		for (int mat1=mat0+1; mat1<meshes.size(); mat1++) {

			cerr<<"triangulating interface: "<<mat0<<" "<<mat1<<endl;
			mmvol->cur_mat1 = mat0;
			mmvol->cur_mat2 = mat1;

			// hack to skip the frog/air interface so we can see the internal components
			// air is given as the 0'th material and only interfaces with the flesh
			// if (mat0==0 || mat1==0)
			// 	continue;

			vector<int> iloops;
			for (int l=0; l<loops.size(); l++) {
				if (contains(loop_mats[l], mat0) &&
					contains(loop_mats[l], mat1)) {
					iloops.push_back(l);
				}
			}

			if (iloops.size() > 0) {
				// we actully have an interface!
				vector< vector<Point3> > loop_pts;
				vector< vector<Vector3> > loop_norms;

				for (int i=0; i<iloops.size(); i++) {
					loop_pts.push_back(vector<Point3>());
					loop_norms.push_back(vector<Vector3>());

					vector<Point3> &loop = loops[iloops[i]];

					for (int j=0; j<loop.size(); j++) {
						Vector3 norm;
						mmvol->NormalAtPoint(loop[j], norm);
						loop_pts.back().push_back(loop[j]);
						loop_norms.back().push_back(norm);
					}

					// if we have just a single point, it's a seed
					if (loop.size() == 1) {

						Point3 p0 = loop[0];
						Vector3 n0 = loop_norms.back()[0];

						Vector3 udir, vdir;
						PerpVectors(n0, udir, vdir);

						Point3 fp = p0 + guidance->MaxStepLength(p0) * udir;

						Point3 tp;
						Vector3 tn;
						projector.ProjectPoint(fp, tp, tn);

						loop_pts.back().push_back(tp);
						loop_norms.back().push_back(tn);

					} else {

						// might need to reverse it so we grow the right way
						// vote for which side we want
						int vote_left = 0;
						int vote_right = 0;

						for (int j=0; j<loop.size(); j++) {
							Vector3 dir = loop_pts.back()[(j+1)%loop_pts.back().size()] - loop_pts.back()[j];
							Vector3 norm = loop_norms.back()[j];

							Vector3 to_right = dir.cross(norm);
							to_right.normalize();
							to_right *= 4*dir.length();

							int mat_right = mmvol->MaterialAtPoint(loop_pts.back()[j] + to_right);
							if (mat_right == mat0 || mat_right == mat1)
								vote_right++;

							int mat_left = mmvol->MaterialAtPoint(loop_pts.back()[j] - to_right);
							if (mat_left == mat0 || mat_left == mat1)
								vote_left++;
						
						}					

						cerr<<vote_left<<endl;
						cerr<<vote_right<<endl;
						if (vote_right > vote_left) {
							reverse(loop_pts.back());
							reverse(loop_norms.back());
						}
					}
				}
				
				
				// run triangulator
				triangulator = new Triangulator(*controller);
				triangulator->SetFlipOutput(true);
				triangulator->SetOutputIndices(num_verts, num_faces);
				triangulator->Go(loop_pts, loop_norms, failsafe);
				triangulator->GetOutputIndices(num_verts, num_faces);
				delete triangulator;
				triangulator = NULL;
			}

		}
	}
	

#if 0
	triangulator = new Triangulator(*controller);
	triangulator->SetFlipOutput(true);

	vector< vector<Point3> > ipts;
	vector< vector<Vector3> > inorms;
	triangulator->Go(ipts, inorms, failsafe);

	cerr<<volume->GetNumBlockLoads()<<" block loads"<<endl;
#endif

	return used_args;

}


int do_multi_material_mc(int argc, char* argv[]) {

	int used_args = 1;
	bool bspline=false;
	if (argc>1 && !stricmp(argv[1], "bspline")) {
		bspline=true;
		used_args++;
	}

	if (triangulator)	delete triangulator;	triangulator=NULL;
	if (guidance)		delete guidance;		guidance=NULL;
	if (controller)		delete controller;		controller=NULL;


	vector<TriangleMesh*> meshes;
	vector<Volume*> volumes;
	for (int i=used_args; i<argc; i+=2) {
		if (argv[i][0] == '-') break;
		meshes.push_back(new TriangleMesh());
		meshes.back()->Read(argv[i]);
		volumes.push_back((Volume*)VolumeRead(argv[i+1]));
		used_args+=2;
	}


	for (int i=0; i<volumes.size(); i++) {

		cerr<<"mc: "<<i<<endl;

		RegularVolume<float> v;
		v.FromMultiMaterial(volumes, i);
		v.SetSplineCoeffs(bspline);

		char name[100];
		sprintf(name, "outmesh.%d.m", i);
		OutputControllerHHM oc(name);

		ControllerWrapper cw(NULL, NULL, &oc);	// for the output
		MarchingCubes(v, cw, -0.0001, true);
	}


	return used_args;
}


class pt_compare {
	public:
	bool operator()(const Point3 &lhs, const Point3 &rhs) {
		if (lhs[0]<rhs[0]) return true;
		if (lhs[0]>rhs[0]) return false;
		if (lhs[1]<rhs[1]) return true;
		if (lhs[1]>rhs[1]) return false;
		if (lhs[2]<rhs[2]) return true;
		if (lhs[2]>rhs[2]) return false;
		return false;
	}
};

void do_draw_material_boundaries() {

	vector< vector<Point3> > loops;


	std::map<Point3, int, pt_compare> ptcount;

	for (int v=0; v<meshes[0].verts.size(); v++) {
		if (ptcount.find(meshes[0].verts[v].point) == ptcount.end()) {
			ptcount[meshes[0].verts[v].point] = 1;
		} else {
			ptcount[meshes[0].verts[v].point]++;
		}
	}

	bool have_interface = false;
	vector<bool> interface(meshes[0].verts.size(), false);

	for (int v=0; v<meshes[0].verts.size(); v++) {
		if (ptcount[meshes[0].verts[v].point] == 2) {
			interface[v] = true;
			have_interface = true;
		}
			
	}


	// pull out the boundary regions of the marked area
	cerr<<"finding boundaries"<<endl;
	bool have_boundary = false;
	vector<bool> done(meshes[0].verts.size(), false);
	for (int v=0; v<meshes[0].verts.size(); v++) {
		if (done[v]) continue;
		if (!interface[v]) {
			done[v] = true;
			continue;
		}

		// see if we're on the boundary
		bool have_marked = false;
		bool have_unmarked = false;
				
		for (TriangleMesh::VertexVertexIteratorI vi(meshes[0], v); !vi.done(); ++vi) {
			int nv = *vi;
			if (interface[nv])
				have_marked = true;
			else
				have_unmarked = true;
		}

		if (!have_marked || !have_unmarked) {
			done[v] = true;
			continue;
		}


		// start a new loop from here
		have_boundary = true;
		//		cerr<<"making loop"<<endl;
		vector<Point3> loop;
		loop.push_back(meshes[0].verts[v].point);
		done[v] = true;

		int lv = v;
				
		while(1) {

			vector<int> ring;
			for (TriangleMesh::VertexVertexIteratorI vi(meshes[0], lv); !vi.done(); ++vi) {
				ring.push_back(*vi);
			}

			int next = -1;
			for (int i=0; i<ring.size(); i++) {
				if (!interface[ring[i]] && interface[ring[(i+1)%ring.size()]]) {
					next = ring[(i+1)%ring.size()];
				}
			}

			if (next == -1) {
				cerr<<"wtf"<<endl;
			}

			if (v == next)
				break;

			lv = next;

			loop.push_back(meshes[0].verts[next].point);
			done[next] = true;

		}


		bool loop_found = false;
		for (int l=0; l<loops.size(); l++) {
			if (contains(loops[l], loop.front())) {
				loop_found = true;
			}
		}

		if (!loop_found) {
			loops.push_back(loop);
		}
	}


	for (int l=0; l<loops.size(); l++) {
		Point3 start = loops[l][0];
		loops[l].push_back(start);
		DbgPLines::add(loops[l], 2, boundary_width);
	}
}


int do_key(int argc, char **argv) {
	if (gui) {
		gui->SendKey(argv[1][0], 0);
	}
	sleep(1);
	return 2;
}


void do_merge_meshes() {

	critical_section->enter();

	int nverts1 = meshes[0].verts.size();
	int nfaces1 = meshes[0].faces.size();

	meshes[0].verts.insert(meshes[0].verts.end(), meshes[1].verts.begin(), meshes[1].verts.end());
	meshes[0].faces.insert(meshes[0].faces.end(), meshes[1].faces.begin(), meshes[1].faces.end());

	for (int v=nverts1; v<meshes[0].verts.size(); v++) {
		if (meshes[0].verts[v].someface >= 0)
			meshes[0].verts[v].someface += nfaces1;
	}

	for (int f=nfaces1; f<meshes[0].faces.size(); f++) {
		for (int i=0; i<3; i++) {
			meshes[0].faces[f].verts[i] += nverts1;

			if (meshes[0].faces[f].nbrs[i] >= 0)
				meshes[0].faces[f].nbrs[i] += nfaces1;
		}
	}


	critical_section->leave();
}


int do_save_ooc(int argc, char* argv[]) {
	char *path = argv[1];
	char *name = argv[2];
	int bsx = atoi(argv[3]);
	int bsy = atoi(argv[4]);
	int bsz = atoi(argv[5]);

	int minx = std::max(atoi(argv[6]),0);
	int miny = std::max(atoi(argv[7]),0);
	int minz = std::max(atoi(argv[8]),0);
	int maxx = std::min(atoi(argv[9]), cvolume->GetDim(0));
	int maxy = std::min(atoi(argv[10]), cvolume->GetDim(1));
	int maxz = std::min(atoi(argv[11]), cvolume->GetDim(2));

	int sampling = atoi(argv[12]);

	if (minx == maxx) {
		minx = 0;
		maxx = cvolume->GetDim(0);
	}
	if (miny == maxy) {
		miny = 0;
		maxy = cvolume->GetDim(1);
	}
	if (minz == maxz) {
		minz = 0;
		maxz = cvolume->GetDim(2);
	}




	if (!cvolume) {
		cerr<<"no volume loaded"<<endl;
		return 7;
	}

	char fullprefix[1024];
	if (endswith(path, "/")) {
		sprintf(fullprefix, "%s%s", path, name);
	} else {
		sprintf(fullprefix, "%s/%s", path, name);
	}

	char nooc_name[1024];
	sprintf(nooc_name, "%s.nooc", fullprefix);

	char blockprefix[1024];
	sprintf(blockprefix, "./%s.blocks/%s", name, name);

	char blockdirectory[1024];
	sprintf(blockdirectory, "%s.blocks/", fullprefix);

#ifdef WIN32
	cerr<<"make sure "<<blockdirectory<<" exists!"<<endl;
#else
	mkdir(blockdirectory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif


	FILE *nooc = fopen(nooc_name, "w");
	if (!nooc) {
		cerr<<"couldn't open "<<nooc_name<<endl;
		return 7;
	}
	fprintf(nooc, "sizes: %d %d %d\n", (maxx-minx)/sampling, (maxy-miny)/sampling, (maxz-minz)/sampling);
	fprintf(nooc, "block size: %d %d %d\n", bsx, bsy, bsz);
	fprintf(nooc, "spacings: %g %g %g\n", cvolume->GetAspect(0)*sampling, cvolume->GetAspect(1)*sampling, cvolume->GetAspect(2)*sampling);
	fprintf(nooc, "block prefix: %s\n", blockprefix);
	fprintf(nooc, "space origin: %g %g %g\n", cvolume->GetTranslation()[0], cvolume->GetTranslation()[1], cvolume->GetTranslation()[2]);
	fclose(nooc);



	for (int bx=0; bx*bsx<(maxx-minx)/sampling; bx++) {
		int xstart = minx/sampling + bx*bsx;
		int xend   = std::min(xstart+bsx, maxx/sampling);
		for (int by=0; by*bsy<(maxy-miny)/sampling; by++) {
			int ystart = miny/sampling + by*bsy;
			int yend   = std::min(ystart+bsy, maxy/sampling);
			for (int bz=0; bz*bsz<(maxz-minz)/sampling; bz++) {
				int zstart = minz/sampling + bz*bsz;
				int zend   = std::min(zstart+bsz, maxz/sampling);

				int block[3] = {bx, by, bz};

				char datablockprefix[1024];
				Volume::BlockString(name, block, datablockprefix);
				(*strrchr(datablockprefix, '.')) = '\0';

				char datafilename[1024];
#ifdef HAS_ZLIB
				sprintf(datafilename, "./%s.gz", datablockprefix);
#else
				sprintf(datafilename, "./%s.raw", datablockprefix);
#endif

				char nhdr_name[1024];
				sprintf(nhdr_name, "%s%s.nhdr", blockdirectory, datablockprefix);

				FILE *nhdr = fopen(nhdr_name, "w");
				if (!nhdr) {
					cerr<<"couldn't open "<<nhdr_name<<endl;
					return 7;
				}

				fprintf(nhdr, "NRRD0001\n");
				fprintf(nhdr, "content: %s\n", name);
				fprintf(nhdr, "dimension: 3\n");
				fprintf(nhdr, "type: %s\n", cvolume->TypeString());
				fprintf(nhdr, "sizes: %d %d %d\n", (xend-xstart), (yend-ystart), (zend-zstart));
				fprintf(nhdr, "spacings: %g %g %g\n", cvolume->GetAspect(0)*sampling, cvolume->GetAspect(1)*sampling, cvolume->GetAspect(2)*sampling);
				fprintf(nhdr, "data file: %s\n", datafilename);
#ifdef HAS_ZLIB
				fprintf(nhdr, "encoding: gzip\n");
#else
				fprintf(nhdr, "encoding: raw\n");
#endif

				Vector3 origin = cvolume->GetTranslation() +
					Vector3(xstart*cvolume->GetAspect(0),
							ystart*cvolume->GetAspect(1),
							zstart*cvolume->GetAspect(2));

				fprintf(nhdr, "space origin: (%g,%g,%g)\n", origin[0], origin[1], origin[2]);
				fclose(nhdr);


				char fulldatafilename[1024];
				sprintf(fulldatafilename, "%s%s", blockdirectory, datafilename);

#ifdef HAS_ZLIB
				gzFile f = gzopen(fulldatafilename, "w");
#else
				FILE *f = fopen(fulldatafilename, "w");
#endif
				if (!f) {
					cerr<<"couldn't open file "<<fulldatafilename<<endl;
				}


				for (int z=zstart; z<zend; z++) {
					for (int y=ystart; y<yend; y++) {
						for (int x=xstart; x<xend; x++) {
							cvolume->WriteValue(x*sampling,y*sampling,z*sampling, f);
						}
					}
				}
		

#ifdef HAS_ZLIB
				gzclose(f);
#else
				fclose(f);
#endif
			}
		}	
	}



	return 13;
}


int do_rho_N(int argc, char* argv[])
{
	if ((argc < 2) || (argv[1][0] == '-'))
		{
			cerr<<"Rho requires a parameter"<<endl;
			return 1;
		}
	rho = (real_type)(M_PI * 2.0 / atof(argv[1]));
	return 2;
}


void do_misc() {
	// miscelaneous function to help debug
}



void do_meshinfo() {

	for (int m=0; m<2; m++) {
		int nbl=-1;	 // num boundary loops
		int nbe=0;	// num boundary edges
		int ncc=-1;	 // num connected components


		if (!cmeshes[m]->verts.size()) continue;

		// check that we're a manifold
		vector<int> vertfaces(cmeshes[m]->verts.size());
		for (unsigned i=0; i<vertfaces.size(); i++) vertfaces[i]=0;

		for (unsigned f=0; f<cmeshes[m]->faces.size(); f++) {
			for (int i=0; i<3; i++) {
				vertfaces[cmeshes[m]->faces[f].verts[i]]++;
			}
		}

		bool manifold=true;
		int numdangling=0;

		// this barfs on non-manifold meshes, so skip it for now!
		manifold=false;
		if (0)
			for (unsigned v=0; v<cmeshes[m]->verts.size(); v++) {

				if (vertfaces[v]==0 || cmeshes[m]->verts[v].someface==-1) {
					cerr<<"dangling vertex (index: "<<v<<")"<<endl;
					numdangling++;
				}

				int tvf = 0;
				for (TriangleMesh::VertexFaceIteratorI vfi(*cmeshes[m], v); !vfi.done(); ++vfi) {
					tvf++;
				}

				if (tvf != vertfaces[v]) {
					cerr<<"non-manifold vertex (index: "<<v<<")"<<endl;
					dbgClear();
					DbgPoints::add(cmeshes[m]->verts[v].point, 1, 0, 0);
					redrawAndWait(' ');
					manifold=false;
				}
			}



		if (!manifold) {
			cerr<<"skipping info for non-manifold surface!"<<endl;
		} else {

			nbl = 0;
			nbe = 0;

			// count the number of boundary loops
			{
				vector<bool> vb(cmeshes[m]->verts.size());
				for (unsigned i=0;i<cmeshes[m]->verts.size(); i++) {
					vb[i] = cmeshes[m]->VertOnBoundary(i);
					if (vb[i]) nbe++;
				}

				while (1) {

					int startv=-1;
					for (unsigned i=0; i<vb.size(); i++) {
						if (vb[i]) {
							startv=i; break;
						}
					}
					if (startv<0) break;

					int bi=startv;
					do {
						vb[bi] = false;
						bi = *(TriangleMesh::VertexVertexIteratorI(*cmeshes[m], bi));
					} while (bi!=startv);

					nbl++;
				}
			}


			// count the number of connected components
			ncc = 0;
			{
				vector<bool> flooded(cmeshes[m]->verts.size());
				for (unsigned i=0;i<cmeshes[m]->verts.size(); i++) {
					flooded[i]=false;
				}


				while (1) {

					int startv=-1;
					for (unsigned i=0; i<flooded.size(); i++) {
						if (!flooded[i]) {
							startv=i; break;
						}
					}
					if (startv<0) break;

					vector<int> toflood;
					toflood.push_back(startv);
					flooded[startv] = true;

					while (!toflood.empty()) {
						int tf = toflood.back();
						toflood.pop_back();
						for (TriangleMesh::VertexVertexIteratorI vi(*cmeshes[m], tf); !vi.done(); ++vi) {
							int tf2 = *vi;
							if (flooded[tf2]) continue;
							flooded[tf2]=true;
							toflood.push_back(tf2);
						}
					}

					ncc++;
				}
		
			}
		}

		int num_edges=0;
		for (unsigned f=0; f<cmeshes[m]->faces.size(); f++) {
			for (int i=0; i<3; i++) {
				if (cmeshes[m]->faces[f].nbrs[i] < (int)f)
					num_edges++;
			}
		}

		int num_verts = cmeshes[m]->verts.size() - numdangling;
		int num_faces = cmeshes[m]->faces.size();

		int genus = (2*ncc - num_verts + num_edges - num_faces - nbl) / 2;


		// compute the smallest angle
		real_type sangle = -1, langle = 1;
		int sav = -1, lav = -1;
		real_type smallest_total_angle = 0.0, largest_total_angle = 0.0;
		for (unsigned f=0; f<cmeshes[m]->faces.size(); f++) {
			float mna = -1.0;
			float mxa = 1.0;
			for (int i=0; i<3; i++) {

				Vector3 e1 = cmeshes[m]->verts[cmeshes[m]->faces[f].verts[(i+1)%3]].point - cmeshes[m]->verts[cmeshes[m]->faces[f].verts[(i+0)%3]].point;
				Vector3 e2 = cmeshes[m]->verts[cmeshes[m]->faces[f].verts[(i+2)%3]].point - cmeshes[m]->verts[cmeshes[m]->faces[f].verts[(i+0)%3]].point;

				real_type ta = e1.dot(e2) / (e1.length() * e2.length());

				if (ta > sangle) {
					sangle = ta;
					sav = cmeshes[m]->faces[f].verts[i];
				}
				if (ta < langle) {
					langle = ta;
					lav = cmeshes[m]->faces[f].verts[i];
				}
				if (ta > mna)
					mna = ta;
				if (ta < mxa)
					mxa = ta;
			}
			smallest_total_angle += mna;
			largest_total_angle += mxa;
		}
		if (sangle < -1) sangle = -1;
		if (sangle > 1)	 sangle = 1;
		if (langle < -1) langle = -1;
		if (langle > 1)	 langle = 1;
		sangle = acos(sangle);
		langle = acos(langle);
		smallest_total_angle = acos(smallest_total_angle / (real_type) num_faces);
		largest_total_angle = acos(largest_total_angle / (real_type) num_faces);
		DbgPoints::add(cmeshes[m]->verts[sav].point, 1,0,0);
		DbgPoints::add(cmeshes[m]->verts[lav].point, 1,0,0);


		cerr<<"info for mesh "<<m<<":"<<endl;
		if (!manifold) cerr<<"NON MANIFOLD MESH"<<endl;
		cerr<<"\tnum vertices: "<<num_verts<<endl;
		cerr<<"\tnum edges: "<<num_edges<<endl;
		cerr<<"\tboundary edges: "<<nbe<<endl;
		cerr<<"\tnum faces: "<<num_faces<<endl;
		cerr<<"\tnum boundary loops: "<<nbl<<endl;
		cerr<<"\tnum connected components: "<<ncc<<endl;
		cerr<<"\tgenus: "<<genus<<endl;
		cerr<<"\tsmallest angle: "<<(sangle * 180/M_PI)<<" degrees"<<endl;
		cerr<<"\tlargest angle: "<<(langle * 180/M_PI)<<" degrees"<<endl;
		cerr<<"\tsmallest angle mean: "<<(smallest_total_angle * (180/M_PI))<<" degrees"<<endl;
		cerr<<"\tlargest angle mean: "<<(largest_total_angle * (180/M_PI))<<" degrees"<<endl;

	}

}

#define NUM_HISTOGRAM_BINS 100
void do_histogram(void) {

	int counts[NUM_HISTOGRAM_BINS];
	for (int i=0; i<NUM_HISTOGRAM_BINS; i++)
		counts[i] = 0;



	if (cmeshes[0]->faces.size()) {
		for (unsigned f=0; f<cmeshes[0]->faces.size(); f++) {

			real_type elen[3] = { Point3::distance(cmeshes[0]->verts[cmeshes[0]->faces[f].verts[0]].point, cmeshes[0]->verts[cmeshes[0]->faces[f].verts[1]].point),
								  Point3::distance(cmeshes[0]->verts[cmeshes[0]->faces[f].verts[1]].point, cmeshes[0]->verts[cmeshes[0]->faces[f].verts[2]].point),
								  Point3::distance(cmeshes[0]->verts[cmeshes[0]->faces[f].verts[2]].point, cmeshes[0]->verts[cmeshes[0]->faces[f].verts[0]].point) };

			real_type circumrad, inrad;
			void circumradius_inradius(real_type a, real_type b, real_type c, real_type &circumrad, real_type &inrad);
			circumradius_inradius(elen[0], elen[1], elen[2], circumrad, inrad);
		
			int bin = (int)(NUM_HISTOGRAM_BINS * ((inrad / circumrad) * 2.0f));
			if (bin < 0) bin = 0;
			if (bin >= NUM_HISTOGRAM_BINS) bin = NUM_HISTOGRAM_BINS-1;
			counts[bin]++;
		}
	} else if (ctetmesh->tets.size()) {
		for (unsigned t=0; t<ctetmesh->tets.size(); t++) {
			real_type ratio = ctetmesh->radius_ratio(t);
			int bin = (int)(NUM_HISTOGRAM_BINS * (ratio * 3.f));
			if (bin < 0) bin = 0;
			if (bin >= NUM_HISTOGRAM_BINS) bin = NUM_HISTOGRAM_BINS-1;
			counts[bin]++;
		}
	}




	int mode = 0;
	for (int i=1; i<NUM_HISTOGRAM_BINS; ++i) {
		if (counts[i] > mode)
			mode = counts[i];
		counts[i] += counts[i-1];
	}
	double quality = 0.0;

	
	if (!histname) {
		return;
	}

	FILE *f = fopen(histname, "w");
	for (int i=0; i<NUM_HISTOGRAM_BINS; i++) {
		fprintf(f, "%g\n", (double)counts[i] / (double) counts[NUM_HISTOGRAM_BINS-1]);
		if (i>0)
			quality += (1.0/99.0) * (counts[i-1]+counts[i])/(2.0 * counts[NUM_HISTOGRAM_BINS-1]);
	}
	fclose(f);
	cerr << "overall mesh quality (1.0 is best): " << 1.0-quality << endl;
		
	if (!histfigname) {
		return;
	}
	FILE *h = fopen(histfigname, "w");
	fprintf(h, "P3 %d %d\n 1\n", NUM_HISTOGRAM_BINS, NUM_HISTOGRAM_BINS);
	for (int i=0; i<NUM_HISTOGRAM_BINS; ++i) {
		// 0<=th<1
		float th = ((double)NUM_HISTOGRAM_BINS - (double)i) / (double) NUM_HISTOGRAM_BINS;
		for (int j=0; j<NUM_HISTOGRAM_BINS; ++j) {
			float count = counts[j]/(float)counts[NUM_HISTOGRAM_BINS-1];
			float singlebin;
			if (j == 0)
				singlebin = counts[0];
			else
				singlebin = counts[j]-counts[j-1];
			singlebin /= (float)mode;

			fprintf(h, "%d 1 %d\n",th>=count ,th>=singlebin);
		}
	}
	fclose(h);
}


void do_dbg() {
	BREAK;
}


void do_sync_cam_out() {
	if (gui) {
		gui->SetCamSyncOut(0);
	}
}

void do_sync_cam_in() {
	if (gui) {
		gui->SetCamSyncIn(0);
	}
}


void read_pc(const char *filename) {
	vector<InfoPoint> pts;
	PC_read(filename, pts);

	for (unsigned i=0; i<pts.size(); i++) {
		points.insert_vertex(Point3(pts[i].x_, pts[i].y_, pts[i].z_));
	}
}


void read_asf(const char *fname, surfel_set &pts) {

	FILE *f = fopen(fname, "r");

	char delims[] = " \t\n{}=()";
	char line[2048];

	while (fgets(line, 2040, f)) {
		Point3 p;
		Vector3 n;
		real_type r;

		p[0] = atof(strtok(line, delims));
		p[1] = atof(strtok(NULL, delims));
		p[2] = atof(strtok(NULL, delims));

		n[0] = atof(strtok(NULL, delims));
		n[1] = atof(strtok(NULL, delims));
		n[2] = atof(strtok(NULL, delims));

		r = atof(strtok(NULL, delims));

		pts.insert_vertex(p);
		pts.insert_normal(n);
		//	pts.insert_radius(r);
	}

	fclose(f);
}


void project_dense_parallel(int nt, int id, SmoothMLSProjector &proj, real_type noise, int totalpts, FILE *f1, FILE *f2) {

	vector<Point3> frompts;
	vector<Point3> fromnorms;
	vector<Point3> topts;
	vector<Vector3> tonorms;

	int cur = 0;

	for (int i=id; i<totalpts; i+=nt) {

		cur += (int)(myran1f(id) * 1024);
		cur %= points.size();

		Point3 from = points.vertex(cur);
		Vector3 fromn = points.normal(cur);
		from += (2*noise) * Vector3(myran1f(id)-0.5,
									myran1f(id)-0.5,
									myran1f(id)-0.5);

		Point3 top;
		Vector3 ton;
	
		if (proj.ProjectPoint(from, top, ton) == PROJECT_SUCCESS) {
			frompts.push_back(from);
			fromnorms.push_back(fromn);
			topts.push_back(top);
			tonorms.push_back(ton);
		}


		if (i/nt % 1000 == 0) {
			static thlib::CSObject cs;
			cs.enter();
			for (unsigned j=0; j<frompts.size(); j++) {
				fprintf(f1, "v %g %g %g\n", topts[j][0], topts[j][1], topts[j][2]);
				fprintf(f1, "vn %g %g %g\n", tonorms[j][0], tonorms[j][1], tonorms[j][2]);

				fprintf(f2, "v %g %g %g\n", frompts[j][0], frompts[j][1], frompts[j][2]);
				fprintf(f2, "vn %g %g %g\n", fromnorms[j][0], fromnorms[j][1], fromnorms[j][2]);
			}
			cs.leave();

			frompts.clear();
			fromnorms.clear();
			topts.clear();
			tonorms.clear();
		}

	}
}

int do_project_dense(int argc, char **argv) {

	real_type noise = atof(argv[1]);
	int totalpts = atoi(argv[2]);

	SmoothMLSProjector mls(points, adamson);
	mls.SetRadiusFactor(radius_factor);
	
	int cur = 0;

	FILE *f = fopen("out.obj", "w");
	FILE *f2 = fopen("from.obj", "w");
	FILE *f3 = fopen("points.tilo", "w");

	/*
	  ParallelExecutor(idealNumThreads, &project_dense_parallel,
	  mls, noise, totalpts, f, f2);

	*/


	while (totalpts-- !=0) {

		cur += (int)(gtb::nrran1f() * 1024);
		cur %= points.size();

		Point3 from = points.vertex(cur);
		from += (2*noise) * Vector3(gtb::nrran1f()-0.5,
									gtb::nrran1f()-0.5,
									gtb::nrran1f()-0.5);

		Point3 top;
		Vector3 ton;
	
		if (mls.ProjectPoint(from, top, ton) == PROJECT_SUCCESS) {
			fprintf(f, "v %g %g %g\n", top[0], top[1], top[2]);
			fprintf(f2, "v %g %g %g\n", from[0], from[1], from[2]);
			fprintf(f3, "%g %g %g %g %g %g 1\n", top[0], top[1], top[2], ton[0], ton[1], ton[2]);
		}
	}


	fclose(f);
	fclose(f2);
	fclose(f3);


	return 3;
}



int do_ppm(int argc, char **argv) {

	char *dir = argv[1];

	char nooc_name[1024];
	if (endswith(dir, "/")) {
		dir[strlen(dir)-1] = 0;
	}

	sprintf(nooc_name, "%s/ppm.nooc", dir);


	FILE *nooc = fopen(nooc_name, "w");
	if (!nooc) {
		cerr<<"couldn't open "<<nooc_name<<endl;
		return 2;
	}
	fprintf(nooc, "sizes: 2048 2048 1920\n");
	fprintf(nooc, "block size: 256 256 128\n");
	fprintf(nooc, "spacings: 1 1 1\n");
	fprintf(nooc, "block prefix: ./ppm_100\n");
	fprintf(nooc, "space origin: 0 0 0\n");
	fclose(nooc);

	int i=0;
	for (int bz=0; bz<15; bz++) {
		for (int by=0; by<8; by++) {
			for (int bx=0; bx<8; bx++) {
		
				int block[3] = { bx, by, bz };

				char oldname[1024];
				sprintf(oldname, "./d_0100_");
				oldname[ 9] = '0' + (i/1000) % 10;
				oldname[10] = '0' + (i/100) % 10;
				oldname[11] = '0' + (i/10) % 10;
				oldname[12] = '0' + (i/1) % 10;
				oldname[13] = 0;

				char newname[1024];
				Volume::BlockString("ppm_100", block, newname);

				char newnamefull[1024];
				sprintf(newnamefull, "%s/%s", dir, newname);


				FILE *nhdr = fopen(newnamefull, "w");
				if (!nhdr) {
					cerr<<"couldn't open "<<newnamefull<<endl;
					return 2;
				}

				fprintf(nhdr, "NRRD0001\n");
				fprintf(nhdr, "content: ppm\n");
				fprintf(nhdr, "dimension: 3\n");
				fprintf(nhdr, "type: uchar\n");
				fprintf(nhdr, "sizes: 256 256 128\n");
				fprintf(nhdr, "spacings: 1 1 1\n");
				fprintf(nhdr, "data file: %s\n", oldname);
				fprintf(nhdr, "encoding: raw\n");

				fclose(nhdr);

				i++;
			}
		}
	}



	return 2;
}

bool read_surface(const char *fname, int i) {
	
	if (endswith(fname, ".m") || endswith(fname, ".off") || endswith(fname, ".sma")) {
		meshes[i].Read(fname);
		meshes[i].BuildStrips(100);
	} else if (endswith(fname, ".vol") ||
			   endswith(fname, ".nhdr") ||
			   endswith(fname, ".nooc")) {
		if (i!=0) {
			cerr<<"can only read one volume!"<<endl;
			return false;
		}
		Volume *rv = VolumeRead(fname);
		volume = rv;
		cvolume= rv;
	} else if (endswith(fname, ".obj")) {
		read_obj(fname, points);
	} else if (endswith(fname, ".asf")) {
		read_asf(fname, points);
	} else if (endswith(fname, ".offt")) {
		tetmesh.ReadOFF(fname);
	} else if (endswith(fname, ".pc")) {
		read_pc(fname);
	} else {
		cerr<<"unknown file type: "<<fname<<endl;
		return false;
	}

	return true;
}



// gcc won't link to these functions if they are not referenced in the main binary
// (they are referenced in gtb)
#ifndef USENR
#include <f2c.h>
extern "C" {
#include <clapack.h>
}
void fix_gcc_linking() {
	ssyev_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	dsyev_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	dgesvd_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
	sgesvd_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}
#endif



int main(int argc, char* argv[]) {


#ifdef WIN32
	SetPriorityClass(GetCurrentProcess(),BELOW_NORMAL_PRIORITY_CLASS);
#endif

	cout.precision(20);
	cerr.precision(20);

	gtb::nrseed(0);


	command_line[0]=0;
	for (int i=1; i<argc; i++) {
		strcat(command_line, argv[i]);
		strcat(command_line, " ");
	}


	gtb::CommandLine cl(argc, argv, "afront <-nogui> mesh1/mls <mesh2> [-commands]");

	CL_ADD_FUN(cl,misc,				  ": miscleaneous");
	CL_ADD_FUN(cl,meshinfo,			  ": display stats about the mesh");
	CL_ADD_FUN(cl,histogram,		  ": write a histogram of the triangle circle ratios to histogram.txt");
	CL_ADD_FUN(cl,quit,				  ": duh");
	CL_ADD_FUN(cl,dbg,				  ": hit a debug break function");
	CL_ADD_FUN(cl,gui,				  ": start the gui");
	CL_ADD_FUN(cl,subdiv,			  "num : apply num iterations of loop subdivision");
	CL_ADD_FUN(cl,trans,			  ": enter mesh transformation mode - for orienting meshes for csg");
	CL_ADD_FUN(cl,flip,				  ": flip the ordering of mesh0 faces");
	CL_ADD_FUN(cl,flip_pts,			  ": flip the orientation of the pointset normals");
	CL_ADD_FUN(cl,compute_pointset_normals,	": compute the normals of the pointset");
	CL_ADD_FUN(cl,flip_components,	  ": allow the user to flip the orientation of each connected componont");
	CL_ADD_FUN(cl,save_r3d,		      "name : write points / normals for input to reconstruct3d");
	CL_ADD_FUN(cl,save_mesh0,		  "name : save mesh[0] to a file");
	CL_ADD_FUN(cl,save_mesh1,		  "name : save mesh[1] to a file");
	CL_ADD_FUN(cl,save_pts,			  "name : save pointset to a file");
	CL_ADD_FUN(cl,save_ooc,			  "path prefix x y z: write out of core volume with block sizes x,y,z");
	CL_ADD_FUN(cl,merge_meshes,		  ": merge mesh1 into mesh0");
	CL_ADD_FUN(cl,multi_material,	  "<bspline> file1.m file1.nhdr ... : triangulate a multimaterial volume set");
	CL_ADD_FUN(cl,multi_material_mc,  "<bspline> file1.m file1.nhdr ... : run mc on a multimaterial volume set");
	CL_ADD_FUN(cl,draw_material_boundaries,  ": draw the boundaries in one of miriah's multi-material meshes");
	CL_ADD_FUN(cl,key,                "k: send key k to the gui");
	CL_ADD_FUN(cl,ppm,				  "path: write the .nooc wrapper around the ppm dataset");
	CL_ADD_VAR(cl,gffile,			  "name : read ideal lengths from this file (only for mls right now)");
	CL_ADD_VAR(cl,rho,				  ": angle subtended on the osculating sphere");
	CL_ADD_FUN(cl,rho_N,			  "N: # of edges that a circle is divided to");
	CL_ADD_VAR(cl,draw_messages,	  ": draw the debug messages on the sceen?	good for screenshots");
	CL_ADD_VAR(cl,boundary_width,	  "w: how thick to draw the boundaries of a mesh");
	CL_ADD_FUN(cl,tet_smooth,		  "N: make N passes of smoothing on the tet mesh scalars");

	CL_ADD_VAR(cl,sharp,			 "ang : cosine of angle between triangle normals for edge to be considered sharp feature (default, -1 - no sharp features are handled)");
	CL_ADD_VAR(cl,small_crease,		  "edge_count : creases with less than <edge_count> edges will be discarded (default 20)");
	CL_ADD_VAR(cl,outname,			  ": the file to output the mesh to");
	CL_ADD_VAR(cl,meshfrontfile,      ": the file name to read the initial front for remeshing from");
	CL_ADD_VAR(cl,noedgeflip,		  ": disable streaming edge flipper on the output stream");
	CL_ADD_VAR(cl,histname,			  ": the file to output the histogram to");
	CL_ADD_VAR(cl,histfigname,		  ": the file to output the histogram figure to");
	CL_ADD_VAR(cl,min_step,			  ": minimum edge length allowed");
	CL_ADD_VAR(cl,max_step,			  ": maximum edge length allowed");
	CL_ADD_VAR(cl,reduction,		  ": percent reduction allowed in a single step");
	CL_ADD_VAR(cl,radius_factor,	  ": set radius factor for mls triangulation");
	CL_ADD_VAR(cl,csg_guidance_blend, ": blend guidance field into existing edge lengths");
	CL_ADD_VAR(cl,failsafe,			  ": close any holes that failed to be triangulated");
	CL_ADD_VAR(cl,feature_interference_hack, ": don't detect interference between \"feature\" fronts - hack for sharp features");
	CL_ADD_VAR(cl,adamson,			  ": use adamson projection instead of standard mls - 0=standard mls, 1=adamson with weighted ave of normals, 2=adamson with covariance normals from weighted average, 3=adamson with covariance normals from point");
	CL_ADD_VAR(cl,nielson,			  ": use nielson tet interpolation");
	CL_ADD_VAR(cl,fence_scale,		  ": scale the height of the fence for front intersection tests");
	CL_ADD_FUN(cl,project_dense,	  "noise num: project a dense set of points onto the mls surface, output in out.obj");

	CL_ADD_FUN(cl,tri_mesh,			  "subdiv : triangulate the mesh, applying subdiv iterations of loop subdivision before computing curvature");
	CL_ADD_FUN(cl,tri_mls,			  ": triangulate a smooth mls surface");
	CL_ADD_FUN(cl,tri_vol,			  "isovalue <bspline>: triangulate an isosurface from the volume");
	CL_ADD_FUN(cl,make_gf,			  "block isovalue <bspline>: create the guidance field for the block num given");
	CL_ADD_FUN(cl,tri_tet,			  "isovalue: triangulate an isosurface from the tet mesh");
	CL_ADD_VAR(cl,stream_mc,          ": use streaming marchingcubes (can't run in parallel)");
	CL_ADD_FUN(cl,marchingcubes,	  "isovalue: extract an isosurface from the regular volume");
	CL_ADD_FUN(cl,volsize,			  ": print out the size of the volume, trying to ignore degenerate isovalues");
	CL_ADD_FUN(cl,segmentation,       ": convert the current volume (segmented) into a volume that represents all the material boundaries");
	CL_ADD_FUN(cl,vol_smooth,         "width: box filter the volume");
	CL_ADD_VAR(cl,mc_stats,			  ": use mc hack for stats");
	CL_ADD_FUN(cl,marchingtets,		  "isovalue: extract an isosurface from the tet volume");
	CL_ADD_VAR(cl,curvature_sub,	  ": how much to divide cells for the guidance field");
	CL_ADD_VAR(cl,eval_sub,			  ": how much to divide cells for evaluating the tetmesh mls isofunction");
	CL_ADD_VAR(cl,iso_scale,		  ": scale volume scalars by this value");
	CL_ADD_VAR(cl,trim_guidance,	  ": trim the guindance field or use the full point set");
	CL_ADD_VAR(cl,trim_bin_size,	  ": when to stop subdivision for guidance field trimming");
    CL_ADD_VAR(cl,gf_prefix,          ": prefix for saving the gf");
    CL_ADD_VAR(cl,reuse_guidance,     ": reuse the already computed guidance field for out of core isosurface extraction");
    CL_ADD_VAR(cl,reuse_seeds,        ": reuse the already computed boundary and seed points for isosurface extraction");
	CL_ADD_VAR(cl,retro,              ": retroactively add points to the guidance field - helps when you dont same the gf densely enough");
	CL_ADD_VAR(cl,value_tol,		  ": the isovalue termination condition for newton stepping onto isosurfaces");
	CL_ADD_VAR(cl,step_tol,			  ": the step size termination condition for newton stepping onto isosurfaces");
	CL_ADD_VAR(cl,max_iter,			  ": maximum number of iterations when projecting points onto isosurfaces");
	CL_ADD_VAR(cl,bad_connect_priority, ": if the triangle priority (circumradius/inradius) is > this, don't consider that triangle");
	CL_ADD_VAR(cl,mls_gf_fix,		  ": iterations for fixing bad curvatures in the gf");
 	CL_ADD_VAR(cl,grow_from_connect,  ": should edge grows be allowed from connection edges");
	CL_ADD_VAR(cl,boundary_dist,	  ": boundary detection parameter for triangle soups - 0 disables (uses topological boundaries), I recommend -0.2.	Negative is relative, positive is absolute");
	CL_ADD_VAR(cl,rmls,				  ": use robust mls instead of standard");
	CL_ADD_VAR(cl,rmls_knn,			  ": number of nearest neighbors to use for rmls projection");
	CL_ADD_VAR(cl,noise_threshold,	  ": noise parameter for rmls projection and feature detection");

	CL_ADD_VAR(cl,max_grab_bag,		   ": the maximum number of blocks allowed in memory for out of core volumes");
	CL_ADD_VAR(cl,prefetch_rings,	   ": number of block rings to prefetch");
	CL_ADD_VAR(cl,prefetch_lookahead,  ": the number of blocks in z-order to look ahead for prefetching");

	CL_ADD_FUN(cl,csg,				  "command <subdiv> : do a csg operation on the two meshes and retriangulate");
	CL_ADD_FUN(cl,tri,				  ": try to smartly decide which triangulation to do");

	CL_ADD_FUN(cl,sync_cam_out,		  "write the camera view to stdout continuously");
	CL_ADD_FUN(cl,sync_cam_in,		  "read the camera view from stdin continuously");

	CL_ADD_VAR(cl,idealNumThreads,	  "num : set the ideal number of execution threads");


	if (argc < 2) {
		cl.Print();
		return 0;
	}

	int toskip=1;

	if ((argc > toskip+1) && (stricmp(argv[toskip], "-nogui")==0))	{
		toskip++;
		critical_section = new thlib::CSObject();
	} else {
		do_gui();
		redrawAll();
		if (gui) gui->ViewAll();
	}


	if ((argc > toskip+1) && (stricmp(argv[toskip], "-noedgeflip")==0))	 {
		noedgeflip = true;
		toskip++;
	}


	critical_section->enter();

	if (argc>toskip && argv[toskip][0] != '-') {
		if (!read_surface(argv[toskip], 0)) return 1;
		toskip++;

		if (argc>toskip && argv[toskip][0] != '-') {
			if (!read_surface(argv[toskip], 1)) return 1;
			toskip++;
		}
	}

	critical_section->leave();
	

	redrawAll();
	if (gui) gui->ViewAll();
	cl.Parse(toskip);

	redrawAndWait('x', true);

	if (gui)
		delete gui;
	return 0;
}


