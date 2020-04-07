/*
===============================================================================

  FILE:  sm_viewer.cpp
  
  CONTENTS:
  
    This little tool visualizes a streaming mesh (and also not so streaming ones).
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    23 May 2007 -- using smreadopener for file IO 
    29 July 2005 -- added support for reading OFF
    24 January 2005 -- improved SME takes place of old SMC
    20 January 2005 -- added out-of-core rendering functionality ('r')
    12 January 2005 -- added support for stdin/stdout
    09 January 2005 -- renamed from 'sm_viz.cpp', added support for SME
    06 January 2004 -- created after 'onyx' burned his tail in my candle
  
===============================================================================
*/
#include <vector>
#include <set>

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>

#include <hash_map.h>

#include <sm/vec3fv.h>
#include <sm/vec3iv.h>

#include <sm/smreadopener.h>

#include <sm/smreadpostascompactpre.h>

#include <viewer/viewer.h>

void MyMenuFunc(int);

void usage()
{
	fprintf(stderr,"usage:\n");
	fprintf(stderr,"sm_viewer mesh.sma\n");
	fprintf(stderr,"sm_viewer -ismc < mesh.smc\n");
	fprintf(stderr,"sm_viewer -win 640 480 \n");
	fprintf(stderr,"sm_viewer -win 1600 1200 mymesh.obj\n");
	fprintf(stderr,"sm_viewer -h\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"in addition to our formats SMA to SMD this tool\n");
	fprintf(stderr,"supports OBJ, SMF, PLY meshes (optionally gzipped).\n");
	exit(0);
}


class SMViewer : public Viewer {

	public:
	SMViewer(int argc, char **argv) {

		// viewer
		char title[1024];
		strcpy(title, "Streaming Mesh Viewer");
		AddWindow(title, 600, 600);



		boundingBoxScale=1.0f;
		boundingBoxTranslateX = 0.0f;
		boundingBoxTranslateY = 0.0f;
		boundingBoxTranslateZ = 0.0f;

		// GLOBAL CONTROL VARIABLES
		doFull = false;
		AnimationOn=0;
		WorkingOn=0;
		PrintToFileOn=0;

		// COLORS
		switch_colors = true;

		// DATA STORAGE FOR STREAM VISUALIZER OUTPUT
		smreadopener = 0;
		smreader = 0;

		max_vertex_span;

		RUN_UNTIL = 0;
		EXACTLY_N_STEPS = 50;
		EVERY_NTH_STEP = -1;
		NEXT_STEP;

		DIRTY_MESH=1;
		DIRTY_VERTEX_SPAN = 1;
		DIRTY_LINES = 1;
		REPLAY_IT=0;
		REPLAY_COUNT=0;
		GRID_PRECISION=3;
		COMPUTE_VERTEX_SPAN=0;
		RED_SCALE = 2;
		STREAM_COLORING = 3;
		NUMBER_STRIPES = 20;
		STRIPE_WHITE=0;
		RENDER_MODE = 0;
		EXTRA_Z_SCALE = 1;

		SHOW_VERTEX_SPAN=0;
		SHOW_LINES=0;

		Framebuffer = 0;
		PrintFileName = "frame";
		Time=0;

		index_map_size = 0;
		index_map_maxsize = 0;

		grid_vertex_buffer = 0;
		grid_vertex_buffer_size = 0;
		grid_vertex_buffer_alloc = 0;

		triangle_buffer = 0;
		triangle_buffer_size = 0;
		triangle_buffer_alloc = 0;

		render_vertex_buffer_alloc = 1024;
		render_vertex_buffer_next = 0;

		

		// first look over the command line arguments for input
		int i;
		for (i = 1; i < argc; i++) {
			if (strncmp(argv[i],"-i",2) == 0) {
				if (smreadopener == 0) {
					// apparently there is a mesh reader to be opened
					smreadopener = new SMreadOpener();
				}
				if (strcmp(argv[i],"-i") == 0) { // the commandline '-i' followed by the input filename
					i++;
					smreadopener->set_file_name(argv[i]);
				}
				else {// the commandline '-ixxx' specifies the input format 'xxx'
					smreadopener->set_file_format(&(argv[i][2]));
				}
			}
		}

		// second look over the command line arguments for options
		for (i = 1; i < argc; i++) {
			if (strncmp(argv[i],"-i",2) == 0) {
				if (strcmp(argv[i],"-i") == 0) i++;
				// already processed
			}
			else if (strcmp(argv[i],"-sco") == 0) {
				i++;
				SetCamSyncOut();
			}
			else if (strcmp(argv[i],"-sci") == 0) {
				i++;
				SetCamSyncIn();
			}
			else if (strcmp(argv[i],"-grid") == 0) {
				i++;
				GRID_PRECISION = atoi(argv[i]);
			}
			else if (strcmp(argv[i],"-every") == 0) {
				i++;
				EVERY_NTH_STEP = atoi(argv[i]);
			}
			else if (strcmp(argv[i],"-until") == 0) {
				i++;
				RUN_UNTIL = atoi(argv[i]);;
			}
			else if (strcmp(argv[i],"-steps") == 0) {
				i++;
				EXACTLY_N_STEPS = atoi(argv[i]);
			}

			else if ((i == argc-1) && smreadopener == 0) {
				// last argument unknown and no input specified yet ... try it
				smreadopener = new SMreadOpener();
				smreadopener->set_file_name(argv[i]);
			}
			else {
				usage();
			}
		}

		// if we had no input try using this as the default input
		if (smreadopener == 0) {
			smreadopener = new SMreadOpener();
			smreadopener->set_file_name("../models/dragon/dragon-spectral.smc");
		}
		
		grid_hash = new my_hash;
		map_hash = new my_hash;
	}


	void Init(int win) {


		// steps sub menu
		int menuSteps = glutCreateMenu(::MyMenuFunc);
		glutAddMenuEntry("in 5 steps", 40);
		glutAddMenuEntry("in 10 steps", 41);
		glutAddMenuEntry("in 25 steps", 42);
		glutAddMenuEntry("in 50 steps", 43);
		glutAddMenuEntry("in 100 steps", 44);
		glutAddMenuEntry("in 250 steps", 45);
		glutAddMenuEntry("in 500 steps", 46);
		glutAddMenuEntry("in 1000 steps", 47);

		// points sub menu
		int menuLines = glutCreateMenu(::MyMenuFunc);
		glutAddMenuEntry("10 lines", 81);
		glutAddMenuEntry("15 lines", 82);
		glutAddMenuEntry("20 lines", 83);
		glutAddMenuEntry("25 lines", 84);
		glutAddMenuEntry("40 lines", 85);
		glutAddMenuEntry("50 lines", 86);
		glutAddMenuEntry("60 lines", 87);
		glutAddMenuEntry("75 lines", 88);

		// points sub menu
		int menuGrid = glutCreateMenu(::MyMenuFunc);
		glutAddMenuEntry("3x3x3 grid", 71);
		glutAddMenuEntry("4x4x4 grid", 72);
		glutAddMenuEntry("5x5x5 grid", 73);
		glutAddMenuEntry("6x6x6 grid", 74);
		glutAddMenuEntry("7x7x7 grid", 75);
		glutAddMenuEntry("8x8x8 grid", 76);
		glutAddMenuEntry("9x9x9 grid", 77);
		glutAddMenuEntry("10x10x10 grid", 78);

		// main menu
		glutCreateMenu(::MyMenuFunc);
		glutAddSubMenu("grid ...", menuGrid);
		glutAddMenuEntry("", 0);
		glutAddSubMenu("steps ...", menuSteps);
		glutAddSubMenu("lines ...", menuLines);
		glutAddMenuEntry("", 0);
		glutAddMenuEntry("<s>tep", 103);
		glutAddMenuEntry("<p>lay", 104);
		glutAddMenuEntry("st<o>p", 105);
		glutAddMenuEntry("", 0);
		glutAddMenuEntry("<v>ertex span (+/-)", 150);
		glutAddMenuEntry("stream-<l>ines", 151);
		glutAddMenuEntry("", 0);
		glutAddMenuEntry("stream <c>oloring", 152);
		glutAddMenuEntry("render <m>ode", 153);
		glutAddMenuEntry("", 0);
		glutAddMenuEntry("<Q>UIT", 109);
		glutAttachMenu(GLUT_RIGHT_BUTTON);



		glShadeModel(GL_FLAT);
  
		InitColors();
		InitLight();

		
		//		find_boundary_verts();

		Key(0,'p'); // play the mesh
		Key(0,'a'); // view all
	}



	void MyMenuFunc(int value)
	{
		if (value >= 100) {
			if (value <= 102) {
				//				InteractionMode = value - 100;
				//				Redraw();
			}
			else if (value == 103) {
				Key(0,'s');
			}
			else if (value == 104) {
				Key(0,'p');
			}
			else if (value == 105) {
				Key(0,'o');
			}
			else if (value == 109) {
				Key(0,'q');
			}
			else if (value == 150) {
				Key(0,'v');
			}
			else if (value == 151) {
				Key(0,'l');
			}
			else if (value == 152) {
				Key(0,'c');
			}
			else if (value == 153) {
				Key(0,'m');
			}
		}

		else if (value == 40) {
			EXACTLY_N_STEPS = 5;
		}
		else if (value == 41) {
			EXACTLY_N_STEPS = 10;
		}
		else if (value == 42) {
			EXACTLY_N_STEPS = 25;
		}
		else if (value == 43) {
			EXACTLY_N_STEPS = 50;
		}
		else if (value == 44) {
			EXACTLY_N_STEPS = 100;
		}
		else if (value == 45) {
			EXACTLY_N_STEPS = 250;
		}
		else if (value == 46) {
			EXACTLY_N_STEPS = 500;
		}
		else if (value == 47) {
			EXACTLY_N_STEPS = 1000;
		}
		else if (value == 71) {
			if (GRID_PRECISION != 3) {
				GRID_PRECISION = 3;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 72) {
			if (GRID_PRECISION != 4) {
				GRID_PRECISION = 4;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 73) {
			if (GRID_PRECISION != 5) {
				GRID_PRECISION = 5;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 74) {
			if (GRID_PRECISION != 6) {
				GRID_PRECISION = 6;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 75) {
			if (GRID_PRECISION != 7) {
				GRID_PRECISION = 7;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 76) {
			if (GRID_PRECISION != 8) {
				GRID_PRECISION = 8;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 77) {
			if (GRID_PRECISION != 9) {
				GRID_PRECISION = 9;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 78) {
			if (GRID_PRECISION != 10) {
				GRID_PRECISION = 10;
				DIRTY_MESH = 1;
				Key(0,'p');
			}
		}
		else if (value == 81) {
			if (NUMBER_STRIPES != 10) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 10;
		}
		else if (value == 82) {
			if (NUMBER_STRIPES != 15) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 15;
		}
		else if (value == 83) {
			if (NUMBER_STRIPES != 20) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 20;
		}
		else if (value == 84) {
			if (NUMBER_STRIPES != 25) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 25;
		}
		else if (value == 85) {
			if (NUMBER_STRIPES != 40) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 40;
		}
		else if (value == 86) {
			if (NUMBER_STRIPES != 50) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 50;
		}
		else if (value == 87) {
			if (NUMBER_STRIPES != 60) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 60;
		}
		else if (value == 88) {
			if (NUMBER_STRIPES != 75) {DIRTY_LINES = 1; SHOW_LINES = 0;}
			NUMBER_STRIPES = 75;
		}

		Redraw();
	}


	// VISUALIZATION SETTINGS
	float boundingBoxScale;
	float boundingBoxTranslateX;
	float boundingBoxTranslateY;
	float boundingBoxTranslateZ;

	// GLOBAL CONTROL VARIABLES
	volatile bool doFull;
	int AnimationOn;
	int WorkingOn;
	int PrintToFileOn;

	// COLORS
	bool switch_colors;
	float colours_diffuse[10][4];
	float colours_white[4];
	float colours_goldenrod[4];
	float colours_light_blue[4];
	float colours_red_white[4];
	float colours_stripe_white[4];
	float colours_stripe_black[4];
	float colours_slow_change[4];

	// DATA STORAGE FOR STREAM VISUALIZER OUTPUT

	SMreadOpener* smreadopener;
	SMreader* smreader;
	double quantizeMultiplier[3];

	int max_vertex_span;

	int RUN_UNTIL;
	int EXACTLY_N_STEPS;
	int EVERY_NTH_STEP;
	int NEXT_STEP;

	int DIRTY_MESH;
	int DIRTY_VERTEX_SPAN;
	int DIRTY_LINES;
	int REPLAY_IT;
	int REPLAY_COUNT;
	int GRID_PRECISION;
	int COMPUTE_VERTEX_SPAN;

	std::vector<float> boundary_verts;
	std::set<int> boundary_vert_indices;

	int RED_SCALE;
	int STREAM_COLORING;
	int NUMBER_STRIPES;
	int STRIPE_WHITE;
	int RENDER_MODE;
	int EXTRA_Z_SCALE;

	int SHOW_VERTEX_SPAN;
	int SHOW_LINES;

	unsigned char* Framebuffer;
	char* PrintFileName;
	int Time;

	// efficient memory allocation

	typedef struct GridVertex
	{
		float v[3];
		int number;
		int index;
		unsigned short v_span;
	} GridVertex;

	typedef struct RenderVertex
	{
		RenderVertex* buffer_next;
		float v[3];
	} RenderVertex;

	typedef struct BoundaryInfo {
		std::vector<int> verts;
		int ntris;
		float pos[3];
	} BoundaryInfo;

	typedef hash_map<int, int> my_hash;
	typedef hash_map<int, RenderVertex*> my_other_hash;
	typedef hash_map<int, BoundaryInfo> my_boundary_hash;

	my_hash* grid_hash;
	my_hash* map_hash;
	int index_map_size;
	int index_map_maxsize;

	GridVertex* grid_vertex_buffer;
	int grid_vertex_buffer_size;
	int grid_vertex_buffer_alloc;


	void initGridVertexBuffer(int alloc)
	{
		if (grid_vertex_buffer) {
			if (grid_vertex_buffer_alloc < alloc) {
				grid_vertex_buffer_alloc = alloc;
				free(grid_vertex_buffer);
				grid_vertex_buffer = (GridVertex*)malloc(sizeof(GridVertex)*grid_vertex_buffer_alloc);
			}
		} else {
			grid_vertex_buffer_alloc = alloc;
			grid_vertex_buffer = (GridVertex*)malloc(sizeof(GridVertex)*grid_vertex_buffer_alloc);
		}
		grid_vertex_buffer_size = 0;
	}


	int allocGridVertex()
	{
		if (grid_vertex_buffer_size == grid_vertex_buffer_alloc) {
			grid_vertex_buffer = (GridVertex*)realloc(grid_vertex_buffer,sizeof(GridVertex)*grid_vertex_buffer_alloc*2);
			if (!grid_vertex_buffer) {
				fprintf(stderr,"FATAL ERROR: realloc grid_vertex_buffer with %d failed.\n",grid_vertex_buffer_alloc*2);
				exit(0);
			}
			grid_vertex_buffer_alloc *= 2;
		}
		int index = grid_vertex_buffer_size;
		grid_vertex_buffer[index].v[0] = 0.0f;
		grid_vertex_buffer[index].v[1] = 0.0f;
		grid_vertex_buffer[index].v[2] = 0.0f;
		grid_vertex_buffer[index].index = index;
		grid_vertex_buffer[index].number = 0;
		grid_vertex_buffer[index].v_span = 0;
		grid_vertex_buffer_size++;
		return index;
	}

	void destroyGridVertexBuffer()
	{
		grid_vertex_buffer_size = 0;
		grid_vertex_buffer_alloc = 0;
		if (grid_vertex_buffer) {
			free(grid_vertex_buffer);
		}
		grid_vertex_buffer = 0;
	}

	int* triangle_buffer;
	int triangle_buffer_size;
	int triangle_buffer_alloc;

	void initTriangleBuffer(int alloc)
	{
		if (triangle_buffer) {
			if (triangle_buffer_alloc < alloc) {
				triangle_buffer_alloc = alloc;
				free(triangle_buffer);
				triangle_buffer = (int*)malloc(sizeof(int)*triangle_buffer_alloc*3);
			}
		} else {
			triangle_buffer_alloc = alloc;
			triangle_buffer = (int*)malloc(sizeof(int)*triangle_buffer_alloc*3);
		}
		triangle_buffer_size = 0;
	}

	int allocTriangle()
	{
		if (triangle_buffer_size == triangle_buffer_alloc) {
			triangle_buffer = (int*)realloc(triangle_buffer,sizeof(int)*triangle_buffer_alloc*3*2);
			if (!triangle_buffer) {
				fprintf(stderr,"FATAL ERROR: realloc triangle_buffer with %d failed.\n",triangle_buffer_alloc*2);
				exit(0);
			}
			triangle_buffer_alloc *= 2;
		}
		int index = triangle_buffer_size;
		triangle_buffer_size++;
		return index;
	}

	void destroyTriangleBuffer()
	{
		triangle_buffer_size = 0;
		triangle_buffer_alloc = 0;
		if (triangle_buffer) {
			free(triangle_buffer);
		}
		triangle_buffer = 0;
	}

	int render_vertex_buffer_alloc;
	RenderVertex* render_vertex_buffer_next;

	RenderVertex* allocRenderVertex(float* v)
	{
		if (render_vertex_buffer_next == 0) {
			render_vertex_buffer_next = (RenderVertex*)malloc(sizeof(RenderVertex)*render_vertex_buffer_alloc);
			if (render_vertex_buffer_next == 0) {
				fprintf(stderr,"malloc for render vertex buffer failed\n");
				return 0;
			}
			for (int i = 0; i < render_vertex_buffer_alloc; i++) {
				render_vertex_buffer_next[i].buffer_next = &(render_vertex_buffer_next[i+1]);
			}
			render_vertex_buffer_next[render_vertex_buffer_alloc-1].buffer_next = 0;
			render_vertex_buffer_alloc = 2*render_vertex_buffer_alloc;
		}
		// get pointer to next available vertex
		RenderVertex* vertex = render_vertex_buffer_next;
		render_vertex_buffer_next = vertex->buffer_next;

		VecCopy3fv(vertex->v, v);
  
		return vertex;
	}

	void deallocRenderVertex(RenderVertex* vertex)
	{
		vertex->buffer_next = render_vertex_buffer_next;
		render_vertex_buffer_next = vertex;
	}

	void SavePPM(char *FileName, unsigned char* Colour, int Width, int Height)
	{
		FILE *fp = fopen(FileName, "wb");
		fprintf(fp, "P6\n%d %d\n255\n", Width, Height);
		int NumRowPixels = Width*3;
		for (int i=(Height-1)*Width*3; i>=0; i-=(Width*3)) {
			fwrite(&(Colour[i]),1,NumRowPixels,fp);
		}
		fclose(fp);
	}

	void InitColors()
	{
		colours_diffuse[0][0] = 0.0f; colours_diffuse[0][1] = 0.0f; colours_diffuse[0][2] = 0.0f; colours_diffuse[0][3] = 1.0f; // black
		colours_diffuse[1][0] = 0.6f; colours_diffuse[1][1] = 0.0f; colours_diffuse[1][2] = 0.0f; colours_diffuse[1][3] = 1.0f; // red
		colours_diffuse[2][0] = 0.0f; colours_diffuse[2][1] = 0.8f; colours_diffuse[2][2] = 0.0f; colours_diffuse[2][3] = 1.0f; // green
		colours_diffuse[3][0] = 0.0f; colours_diffuse[3][1] = 0.0f; colours_diffuse[3][2] = 0.6f; colours_diffuse[3][3] = 1.0f; // blue
		colours_diffuse[4][0] = 0.6f; colours_diffuse[4][1] = 0.6f; colours_diffuse[4][2] = 0.0f; colours_diffuse[4][3] = 1.0f; // yellow
		colours_diffuse[5][0] = 0.6f; colours_diffuse[5][1] = 0.0f; colours_diffuse[5][2] = 0.6f; colours_diffuse[5][3] = 1.0f; // purple
		colours_diffuse[6][0] = 0.0f; colours_diffuse[6][1] = 0.6f; colours_diffuse[6][2] = 0.6f; colours_diffuse[6][3] = 1.0f; // cyan
		colours_diffuse[7][0] = 0.7f; colours_diffuse[7][1] = 0.7f; colours_diffuse[7][2] = 0.7f; colours_diffuse[7][3] = 1.0f; // white
		colours_diffuse[8][0] = 0.2f; colours_diffuse[8][1] = 0.2f; colours_diffuse[8][2] = 0.6f; colours_diffuse[8][3] = 1.0f; // light blue
		colours_diffuse[9][0] = 0.9f; colours_diffuse[9][1] = 0.4f; colours_diffuse[9][2] = 0.7f; colours_diffuse[9][3] = 1.0f; // violett
  
		colours_white[0] = 0.7f; colours_white[1] = 0.7f; colours_white[2] = 0.7f; colours_white[3] = 1.0f; // white
		colours_goldenrod[0] = 205.f/255.f; colours_goldenrod[1] = 190.f/255.f; colours_goldenrod[2] = 112.f/255.f; colours_goldenrod[3] = 1.0f; // lightgoldenrod3
		colours_light_blue[0] = 0.2f; colours_light_blue[1] = 0.2f; colours_light_blue[2] = 0.6f; colours_light_blue[3] = 1.0f; // light blue
		colours_red_white[0] = 0.7f; colours_red_white[1] = 0.7f; colours_red_white[2] = 0.7f; colours_red_white[3] = 1.0f; // white or red

		colours_stripe_white[0] = 0.8f; colours_stripe_white[1] = 0.8f; colours_stripe_white[2] = 0.8f; colours_stripe_white[3] = 1.0f; // white
		colours_stripe_black[0] = 0.5f; colours_stripe_black[1] = 0.5f; colours_stripe_black[2] = 0.5f; colours_stripe_black[3] = 1.0f; // black

		colours_slow_change[0] = 0.5f; colours_slow_change[1] = 0.5f; colours_slow_change[2] = 0.5f; colours_stripe_black[3] = 1.0f; // black
	} 

	void InitLight()
	{
		float intensity[] = {1,1,1,1};
		float position[] = {1,1,5,0}; // directional behind the viewer
		glLightfv(GL_LIGHT0,GL_DIFFUSE,intensity);
		glLightfv(GL_LIGHT0,GL_SPECULAR,intensity);
		glLightfv(GL_LIGHT0,GL_POSITION,position);
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
	}

	void InitGrid()
	{
		for (int i = 0; i < 3; i++) {
			double range = (double)(smreader->bb_max_f[i])-(double)(smreader->bb_min_f[i]);
			if (range > 0.0) {
				quantizeMultiplier[i] = (double)((1<<GRID_PRECISION)-1)/range;
			} else {
				quantizeMultiplier[i] = 0.0;
			}
		}
	}

	int GetGridIndex(float* position)
	{
		int grid_pos[3];
		for (int i = 0; i < 3; i++) {
			grid_pos[i] = (int)(quantizeMultiplier[i] * ((double)position[i] - (double)(smreader->bb_min_f[i])) + 0.5); 
		}
		return (grid_pos[0] << (GRID_PRECISION+GRID_PRECISION)) + (grid_pos[1] << GRID_PRECISION) + (grid_pos[2]);
	}


	void vizBegin()
	{
		REPLAY_IT = 0; // just making sure
		DIRTY_MESH = 1;
		DIRTY_LINES = 1;
		DIRTY_VERTEX_SPAN = 1;

		max_vertex_span = 0;

		if (smreadopener->file_name) {
			fprintf(stderr,"loading '%s' ...\n", smreadopener->file_name);
		}

		if (smreader) { // if we have a mesh reader then just re-open the file
			if (!(smreadopener->reopen(smreader))) {
				fprintf(stderr,"ERROR: cannot reopen smreader\n");
				usage();
			}
		} else {
			if (!(smreader = smreadopener->open())) {
				fprintf(stderr,"ERROR: cannot open the mesh reader\n");
				usage();
			}
			if (smreader->bb_min_f == 0 || smreader->bb_max_f == 0) {
				fprintf(stderr,"need additional pass to compute bounding box ...\n");
				if (smreader->compute_bounding_box()) {
					if (!(smreadopener->reopen(smreader))) {
						fprintf(stderr,"ERROR: cannot reopen after computing the bounding box\n");
						usage();
					}
				}
			}
		}

		if (smreader->post_order) {
			fprintf(stderr,"WARNING: input streaming mesh is post-order\n");
			fprintf(stderr,"         using 'readPostAsCompactPre' ...\n");
			SMreadPostAsCompactPre* smreadpostascompactpre = new SMreadPostAsCompactPre();
			smreadpostascompactpre->open(smreader);
			smreader = smreadpostascompactpre;
		}



		InitGrid();

		fprintf(stderr, "bounding box: x[%g, %g] y[%g, %g] z[%g, %g]", smreader->bb_min_f[0], smreader->bb_max_f[0], smreader->bb_min_f[1], smreader->bb_max_f[1], smreader->bb_min_f[2], smreader->bb_max_f[2]);

		if ((smreader->bb_max_f[1]-smreader->bb_min_f[1]) > (smreader->bb_max_f[0]-smreader->bb_min_f[0])) {
			if ((smreader->bb_max_f[1]-smreader->bb_min_f[1]) > (smreader->bb_max_f[2]-smreader->bb_min_f[2])) {
				boundingBoxScale = 1.0f/(smreader->bb_max_f[1]-smreader->bb_min_f[1]);
			} else {
				boundingBoxScale = 1.0f/(smreader->bb_max_f[2]-smreader->bb_min_f[2]);
			}
		} else {
			if ((smreader->bb_max_f[0]-smreader->bb_min_f[0]) > (smreader->bb_max_f[2]-smreader->bb_min_f[2])) {
				boundingBoxScale = 1.0f/(smreader->bb_max_f[0]-smreader->bb_min_f[0]);
			} else {
				boundingBoxScale = 1.0f/(smreader->bb_max_f[2]-smreader->bb_min_f[2]);
			}
		}
		boundingBoxTranslateX = - boundingBoxScale * (smreader->bb_min_f[0] + 0.5f * (smreader->bb_max_f[0]-smreader->bb_min_f[0]));
		boundingBoxTranslateY = - boundingBoxScale * (smreader->bb_min_f[1] + 0.5f * (smreader->bb_max_f[1]-smreader->bb_min_f[1]));
		boundingBoxTranslateZ = - boundingBoxScale * (smreader->bb_min_f[2] + 0.5f * (smreader->bb_max_f[2]-smreader->bb_min_f[2]));

		fprintf(stderr,"nverts %d nfaces %d\n", smreader->nverts, smreader->nfaces);

		if (EVERY_NTH_STEP == -1) {
			if (smreader->nfaces == -1) {
				EVERY_NTH_STEP = 10000;
			} else {
				EVERY_NTH_STEP = smreader->nfaces / EXACTLY_N_STEPS;
			}
		}
		if (EVERY_NTH_STEP == 0) {
			EVERY_NTH_STEP = 1;
		}
		NEXT_STEP = EVERY_NTH_STEP;

		if (RUN_UNTIL) {
			NEXT_STEP = RUN_UNTIL;
		}

		initGridVertexBuffer(65536);
		initTriangleBuffer(65536*2);

		grid_hash->clear();
		map_hash->clear();
		index_map_maxsize = 0;
		index_map_size = 0;
	}

	void vizEnd()
	{
		REPLAY_IT = 0; // just making sure
		REPLAY_COUNT = triangle_buffer_size;
		DIRTY_MESH = 0;
		DIRTY_LINES = 1;
		if (COMPUTE_VERTEX_SPAN) {
			DIRTY_VERTEX_SPAN = 0;
			COMPUTE_VERTEX_SPAN = 0;
			fprintf(stderr,"... vertex_span is %2.2f %%\n", 100.0f*max_vertex_span/65536);
		}
			
		fprintf(stderr,"grid: %dx%dx%d\n", GRID_PRECISION, GRID_PRECISION, GRID_PRECISION);
		fprintf(stderr,"  simplified vertices: %d\n", grid_vertex_buffer_size);
		fprintf(stderr,"  simplified triangles: %d\n", triangle_buffer_size);
		fprintf(stderr,"maxsize of index hash_map: %d (i.e. front width %2.2f %%)\n", index_map_maxsize, ((float)index_map_maxsize)*100/smreader->v_count);
		index_map_maxsize = 0;
		if (map_hash->size()) {
			fprintf(stderr,"EFFICIENCY WARNING: %d unfinalized streaming mesh vertices\n",map_hash->size());
			map_hash->clear();
		}
	}

	int vizContinue()
	{
		SMevent element = SM_ERROR;
		my_hash::iterator hash_element;

		REPLAY_IT = 0; // just making sure

		while (element) {
			element = smreader->read_element();

			if (element == SM_TRIANGLE) {
				int i;
				int t_idx[3];
				for (i = 0; i < 3; i++) {
					hash_element = map_hash->find(smreader->t_idx[i]);
					if (hash_element == map_hash->end()) {
						fprintf(stderr,"FATAL ERROR: index not in index map hash.\n");
						fprintf(stderr,"             either the streaming mesh is not in\n");
						fprintf(stderr,"             pre-order or it is corrupt.\n");
						exit(0);
					} else {
						t_idx[i] = (*hash_element).second;
					}
					if (smreader->t_final[i]) {
						map_hash->erase(hash_element);
						index_map_size--;
					}
				}
				if (COMPUTE_VERTEX_SPAN) {
					unsigned short v_span = (unsigned short)((((float)(VecMax3iv(smreader->t_idx) - VecMin3iv(smreader->t_idx)))/((float)smreader->nverts))*65535);
					if (v_span > max_vertex_span) max_vertex_span = v_span;

					for (i = 0; i < 3; i++) {
						if (v_span > (grid_vertex_buffer[t_idx[i]].v_span)) {
							grid_vertex_buffer[t_idx[i]].v_span = v_span;
						}
					}
				}
				if (t_idx[0] != t_idx[1] && t_idx[0] != t_idx[2] && t_idx[1] != t_idx[2]) {
					int triangle_idx = allocTriangle();
					VecCopy3iv(&(triangle_buffer[triangle_idx*3]),t_idx);
				}
			} else if (element == SM_VERTEX) {
				int grid_idx, vertex_idx;

				grid_idx = GetGridIndex(smreader->v_pos_f);

				// check if grid vertex for the grid cell of this vertex already exists
				hash_element = grid_hash->find(grid_idx);
				if (hash_element == grid_hash->end()) {
					vertex_idx = allocGridVertex();
					// all following vertices falling into this grid cell find their grid vertex in the grid_hash
					grid_hash->insert(my_hash::value_type(grid_idx, vertex_idx));
				} else {
					vertex_idx = (*hash_element).second;
				}

				// all following triangles can find this vertex in the map_hash
				map_hash->insert(my_hash::value_type(smreader->v_idx, vertex_idx));
				index_map_size++;
				if (index_map_size > index_map_maxsize) index_map_maxsize = index_map_size;

				VecSelfAdd3fv(grid_vertex_buffer[vertex_idx].v, smreader->v_pos_f);
				grid_vertex_buffer[vertex_idx].number++;
			} else if (element == SM_FINALIZED) {
				hash_element = map_hash->find(smreader->final_idx);
				if (hash_element == map_hash->end()) {
					fprintf(stderr,"FATAL ERROR: non-existing vertex was explicitely finalized.\n");
					fprintf(stderr,"             the streaming mesh is corrupt.\n");
					exit(0);
				} else {
					map_hash->erase(hash_element);
					index_map_size--;
				}
			} else if (element == SM_EOF) {
				index_map_size = 0;
				grid_hash->clear();
				map_hash->clear();
			}

			if (smreader->f_count > NEXT_STEP) {
				NEXT_STEP += EVERY_NTH_STEP;
				break;
			}
		}
		if (element) {
			return 1;
		} else {
			return 0;
		}
	}


	// void myIdle()
	// {
	// 	if (AnimationOn) {
	// 		AnimationOn = vizContinue();
	// 		if (!AnimationOn) {
	// 			WorkingOn = 0;
	// 			vizEnd();
	// 		}
	// 		glutPostRedisplay();
	// 	} else if (REPLAY_IT) {
	// 		REPLAY_COUNT += NEXT_STEP;
	// 		glutPostRedisplay();
	// 	}
	// }

	void compute_lines()
	{
		int i;

		if (triangle_buffer_size) {
			for (i = 0; i < grid_vertex_buffer_size; i++) {
				grid_vertex_buffer[i].v_span = 0;
			}
			for (i = 0; i < triangle_buffer_size; i++) {
				if ((i % (triangle_buffer_size/NUMBER_STRIPES)) == 0) {
					STRIPE_WHITE = !STRIPE_WHITE;
				}
				if (STRIPE_WHITE) {
					for (int j = 0; j < 3; j++) grid_vertex_buffer[triangle_buffer[3*i+j]].v_span |= 1;
				} else {
					for (int j = 0; j < 3; j++) grid_vertex_buffer[triangle_buffer[3*i+j]].v_span |= 2;
				}
			}
		}
	}

	void full_resolution_rendering()
	{
		int f_count = 2000000000;

		if (smreadopener->file_name == 0) {
			fprintf(stderr,"ERROR: no input file\n");
			return;
		}

		if (smreader) {
			f_count = smreader->f_count;
			fprintf(stderr,"out-of-core rendering of %d mesh faces ... \n",f_count);
			if (!(smreadopener->reopen(smreader))) exit(1);
		} else {
			fprintf(stderr,"out-of-core rendering of mesh ... \n");
			if (!(smreader = smreadopener->open())) exit(1);
		}

		my_other_hash* render_hash = new my_other_hash;

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
		// glViewport(0,0,WindowW,WindowH);

		// glMatrixMode(GL_PROJECTION);
		// glLoadIdentity();
		// gluPerspective(30.0f,(float)WindowW/WindowH,0.0625f,5.0f);

		glMatrixMode(GL_MODELVIEW);
		// glLoadIdentity();
		// gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0);

		// glRotatef(Elevation,1,0,0);
		// glRotatef(Azimuth,0,1,0);

		// glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ*EXTRA_Z_SCALE);
		// glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale*EXTRA_Z_SCALE);


		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_NORMALIZE);

		if (RENDER_MODE == 0 || RENDER_MODE==3) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else if (RENDER_MODE == 1) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else if (RENDER_MODE == 2) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		}

		if (switch_colors) {
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_white);
			//			glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_goldenrod);
			glMaterialfv(GL_BACK, GL_DIFFUSE, colours_light_blue);
		}
		else {
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_light_blue);
			glMaterialfv(GL_BACK, GL_DIFFUSE, colours_white);
			//			glMaterialfv(GL_BACK, GL_DIFFUSE, colours_goldenrod);
		}

		if (f_count == 0) {
			f_count = smreader->nfaces;
		}

		if (RENDER_MODE == 0) {
			glDisable(GL_COLOR_MATERIAL);
			glShadeModel(GL_FLAT);
		} else {
			glEnable(GL_COLOR_MATERIAL);
			glShadeModel(GL_SMOOTH);
		}

		int i;
		float n[3];
		RenderVertex* v[3];
		my_other_hash::iterator hash_element;
		glBegin(GL_TRIANGLES);
		while (smreader->f_count < f_count) {
			switch (smreader->read_element()) {
			case SM_TRIANGLE:
				for (i = 0; i < 3; i++) {
					hash_element = render_hash->find(smreader->t_idx[i]);
					if (hash_element == render_hash->end()) {
						fprintf(stderr,"FATAL ERROR: index not in render map hash.\n");
						fprintf(stderr,"             either the streaming mesh is not in\n");
						fprintf(stderr,"             pre-order or it is corrupt.\n");
						exit(0);
					}
					else {
						v[i] = (*hash_element).second;
					}
					if (smreader->t_final[i]) {
						deallocRenderVertex(v[i]);
						render_hash->erase(hash_element);
					}
				}

				if (RENDER_MODE == 0) {
					VecCcwNormal3fv(n, v[0]->v, v[1]->v, v[2]->v);
					glNormal3fv(n);
					glVertex3fv(v[0]->v);
					glVertex3fv(v[1]->v);
					glVertex3fv(v[2]->v);
				} else if (boundary_vert_indices.find(smreader->t_idx[0]) != boundary_vert_indices.end() ||
						   boundary_vert_indices.find(smreader->t_idx[1]) != boundary_vert_indices.end() ||
						   boundary_vert_indices.find(smreader->t_idx[2]) != boundary_vert_indices.end()) {

					VecCcwNormal3fv(n, v[0]->v, v[1]->v, v[2]->v);
					glNormal3fv(n);

					for (int i=0; i<3; i++) {
						if (boundary_vert_indices.find(smreader->t_idx[i]) != boundary_vert_indices.end()) {
							glColor3f(1,0.5,0.5);
						} else {
							glColor3f(0.5,0.5,1);
						}
						glVertex3fv(v[i]->v);
					}
				}

				break;
			case SM_VERTEX:
				v[0] = allocRenderVertex(smreader->v_pos_f);
				render_hash->insert(my_other_hash::value_type(smreader->v_idx, v[0]));
				break;
			case SM_FINALIZED:
				hash_element = render_hash->find(smreader->final_idx);
				if (hash_element == render_hash->end()) {
					fprintf(stderr,"FATAL ERROR: index not in render map hash.\n");
					fprintf(stderr,"             either the streaming mesh is not in\n");
					fprintf(stderr,"             pre-order or it is corrupt.\n");
					exit(0);
				}
				deallocRenderVertex((*hash_element).second);
				render_hash->erase(hash_element);
				break;
			case SM_EOF:
				f_count = 0;
				break;
			}
		}
		glEnd();

		// delete any remaining vertices
		while (render_hash->size()) {
			hash_element = render_hash->begin();
			deallocRenderVertex((*hash_element).second);
			render_hash->erase(hash_element);
		}
		delete render_hash;

		fprintf(stderr,"done\n");


		SetUnsharpWidth(25);
		SetUnsharpLambda(30);
		ApplyUnsharp();
		SetUnsharpWidth(0);
		SetUnsharpLambda(4);


		ScreenShot("fullres.ppm", false);
	}


	void find_boundary_verts()
	{
		int f_count = 2000000000;


		if (smreadopener->file_name) {
			fprintf(stderr,"loading '%s' ...\n", smreadopener->file_name);
		}

		if (smreader) { // if we have a mesh reader then just re-open the file
			if (!(smreadopener->reopen(smreader))) {
				fprintf(stderr,"ERROR: cannot reopen smreader\n");
				usage();
			}
		} else {
			if (!(smreader = smreadopener->open())) {
				fprintf(stderr,"ERROR: cannot open the mesh reader\n");
				usage();
			}
			if (smreader->bb_min_f == 0 || smreader->bb_max_f == 0) {
				fprintf(stderr,"need additional pass to compute bounding box ...\n");
				if (smreader->compute_bounding_box()) {
					if (!(smreadopener->reopen(smreader))) {
						fprintf(stderr,"ERROR: cannot reopen after computing the bounding box\n");
						usage();
					}
				}
			}
		}

		if (smreader->post_order) {
			fprintf(stderr,"WARNING: input streaming mesh is post-order\n");
			fprintf(stderr,"         using 'readPostAsCompactPre' ...\n");
			SMreadPostAsCompactPre* smreadpostascompactpre = new SMreadPostAsCompactPre();
			smreadpostascompactpre->open(smreader);
			smreader = smreadpostascompactpre;
		}



		// if (smreadopener->file_name == 0) {
		// 	fprintf(stderr,"ERROR: no input file\n");
		// 	return;
		// }

		// if (smreader) {
		// 	f_count = smreader->f_count;
		// 	fprintf(stderr,"finding boundary of %d mesh faces ... \n",f_count);
		// 	if (!(smreadopener->reopen(smreader))) exit(1);
		// } else {
		// 	fprintf(stderr,"finding boundary of mesh ... \n");
		// 	if (!(smreader = smreadopener->open())) exit(1);
		// }

		// FILE *bmesh = fopen("bmesh.m", "w");
		// int bmesh_v = 1;
		// int bmesh_t = 1;

		my_boundary_hash* render_hash = new my_boundary_hash;

		if (f_count == 0) {
			f_count = smreader->nfaces;
		}

		int i;
		my_boundary_hash::iterator hash_element;
		while (smreader->f_count < f_count) {
			switch (smreader->read_element()) {
			case SM_TRIANGLE:

				for (i = 0; i < 3; i++) {
					hash_element = render_hash->find(smreader->t_idx[i]);
					if (hash_element == render_hash->end()) {
						fprintf(stderr,"FATAL ERROR: index not in render map hash.\n");
						fprintf(stderr,"             either the streaming mesh is not in\n");
						fprintf(stderr,"             pre-order or it is corrupt.\n");
						exit(0);
					}

					BoundaryInfo &bi = (*hash_element).second;
					bi.ntris++;
					for (int j=1; j<3; j++) {
						int vi = smreader->t_idx[(i+j)%3];

						bool found = false;
						for (int k=0; k<bi.verts.size(); k++) {
							if (bi.verts[k] == vi) {
								found=true;
								break;
							}
						}
						if (!found) {
							bi.verts.push_back(vi);
						}
					}

					if (smreader->t_final[i]) {
						//						fprintf(stderr, "final-t %d\n", smreader->t_idx[i]);
						if (bi.verts.size() != bi.ntris) {
							for (int j=0; j<3; j++) {
								boundary_verts.push_back(bi.pos[j]);
							}
							boundary_verts.push_back(bi.ntris);
							boundary_vert_indices.insert(smreader->t_idx[i]);
						}
						render_hash->erase(hash_element);
					}

					// fprintf(bmesh, "Vertex  %d  %g %g %g\n", bmesh_v+i, bi.pos[0], bi.pos[1], bi.pos[2]);
				}

				// fprintf(bmesh, "Face  %d  %d %d %d\n", bmesh_t, bmesh_v+0, bmesh_v+1, bmesh_v+2);

				// bmesh_v += 3;
				// bmesh_t += 1;

				break;

			case SM_VERTEX:
				{
					//					fprintf(stderr, "new vert: %d\n", smreader->v_idx);
					BoundaryInfo bi;
					bi.ntris=0;
					bi.pos[0] = smreader->v_pos_f[0];
					bi.pos[1] = smreader->v_pos_f[1];
					bi.pos[2] = smreader->v_pos_f[2];
					render_hash->insert(my_boundary_hash::value_type(smreader->v_idx, bi));
				}
				break;

			case SM_FINALIZED:
				{
					//					fprintf(stderr, "final-v %d\n", smreader->final_idx);
					hash_element = render_hash->find(smreader->final_idx);
					if (hash_element == render_hash->end()) {
						fprintf(stderr,"FATAL ERROR: index not in render map hash.\n");
						fprintf(stderr,"             either the streaming mesh is not in\n");
						fprintf(stderr,"             pre-order or it is corrupt.\n");
						exit(0);
					}

					BoundaryInfo &bi = (*hash_element).second;
					if (bi.verts.size() != bi.ntris) {
						for (int j=0; j<3; j++) {
							boundary_verts.push_back(bi.pos[j]);
						}
						boundary_verts.push_back(bi.ntris);
						boundary_vert_indices.insert(smreader->final_idx);
					}
					render_hash->erase(hash_element);
				}
				break;

			case SM_EOF:
				f_count = 0;
				break;
			}
		}

		// delete any remaining vertices
		while (render_hash->size()) {
			hash_element = render_hash->begin();

			BoundaryInfo &bi = (*hash_element).second;
			if (bi.verts.size() != bi.ntris) {
				//				fprintf(stderr, "vert %d\n", (*hash_element).first);
				//				fprintf(stderr, "nverts: %d, ntris: %d\n", bi.verts.size(), bi.ntris);
				for (int j=0; j<3; j++) {
					boundary_verts.push_back(bi.pos[j]);
				}
				boundary_verts.push_back(bi.ntris);
				boundary_vert_indices.insert((*hash_element).first);
			}
			render_hash->erase(hash_element);
		}
		delete render_hash;

		// fclose(bmesh);
		fprintf(stderr,"done\n");

		// FILE *pts = fopen("bpts.obj", "w");
		// for (int v=0; v<boundary_verts.size(); v+=4) {
		// 	fprintf(pts, "v %g %g %g\n", boundary_verts[v], boundary_verts[v+1], boundary_verts[v+2]);
		// }		
		// fclose(pts);

	}



	void Mouse(int win, int button, int motion, int wx, int wy, float depth, float x, float y, float z) {
		if (button == GLUT_MIDDLE_BUTTON && motion==-1) {
			if (z < 1e8) {
				LookAt(win, x, y, z);
				cerr<<"looking at "<<x<<" "<<y<<" "<<z<<endl;
				cerr<<"looking at "<<x/256<<" "<<y/256<<" "<<z/128<<endl;
			} else {
				cerr<<"Center not set"<<endl;
			}
		}
		Redraw();
	}
	void PassiveMotion(int win, int wx, int wy, float x, float y, float z) {
	}


	void Idle() {
		if (AnimationOn) {
			AnimationOn = vizContinue();
			if (!AnimationOn) {
				WorkingOn = 0;
				vizEnd();
			}
			Redraw();
		}
		else if (REPLAY_IT) {
			REPLAY_COUNT += NEXT_STEP;
			Redraw();
		}
	}

	void Key(int win, unsigned char key) {

		switch(key)
			{
			case 'a':
			case 'A':
				ViewAll();
				break;
			case 'Q':
			case 'q':
			case 27:
				exit(0);
				break;
			case '-':
				RED_SCALE /= 2;
				if (RED_SCALE < 1) RED_SCALE = 1;
				fprintf(stderr,"RED_SCALE %d\n",RED_SCALE);
				break;
			case '=':
				RED_SCALE *= 2;
				fprintf(stderr,"RED_SCALE %d\n",RED_SCALE);
				break;
			case '[':
				if (EXTRA_Z_SCALE > 1)
					{
						EXTRA_Z_SCALE = EXTRA_Z_SCALE >> 1;
					}
				else
					{
						EXTRA_Z_SCALE = 0;
					}
				fprintf(stderr,"EXTRA_Z_SCALE %d\n",EXTRA_Z_SCALE);
				break;
			case ']':
				if (EXTRA_Z_SCALE > 0)
					{
						EXTRA_Z_SCALE = EXTRA_Z_SCALE << 1;
					}
				else
					{
						EXTRA_Z_SCALE = 1;
					}
				fprintf(stderr,"EXTRA_Z_SCALE %d\n",EXTRA_Z_SCALE);
				break;
			case 'C':
			case 'c':
				STREAM_COLORING = (STREAM_COLORING+1)%4;
				fprintf(stderr,"STREAM_COLORING %d\n",STREAM_COLORING);
				break;
			case 'M':
			case 'm':
				RENDER_MODE = (RENDER_MODE+1)%4;
				fprintf(stderr,"RENDER_MODE %d\n",RENDER_MODE);
				break;
			case 'O':
			case 'o':
				AnimationOn = 0;
				REPLAY_IT = 0;
				break;
			case 'V':
			case 'v':
				fprintf(stderr,"SHOW_VERTEX_SPAN %d\n",!SHOW_VERTEX_SPAN);
				if (SHOW_VERTEX_SPAN)
					{
						SHOW_VERTEX_SPAN = 0;
					}
				else
					{
						SHOW_VERTEX_SPAN = 1;
						if (DIRTY_VERTEX_SPAN)
							{
								COMPUTE_VERTEX_SPAN = 1;
								DIRTY_MESH = 1;
								fprintf(stderr,"computing vertex_span ...\n");
								Key(win, 'p');
							}
					}
				break;
			case 'R':
			case 'r':
				doFull=true;
				break;
			case 'W':
				break;
			case 'w':
				switch_colors = !switch_colors;
				break;
			case 'D':
			case 'd':
				//PrintToFileOn=1;
				//if (Framebuffer) delete [] Framebuffer;
				//Framebuffer = new unsigned char[WindowW*WindowH*3];
				fprintf(stderr,"print_to_file %d\n",PrintToFileOn);
				break;
			case 'T':
				if (DIRTY_MESH)
					{
						// works only in replay mode
						fprintf(stderr,"tiny steps only work during second play (replay)\n");
					}
				else
					{
						REPLAY_COUNT -= 1;
						if (REPLAY_COUNT < 0)
							{
								REPLAY_COUNT = 0;
							}
					}
				break;
			case 't':
				if (DIRTY_MESH)
					{
						// works only in replay mode
						fprintf(stderr,"tiny steps only work during second play (replay)\n");
					}
				else
					{
						if (REPLAY_COUNT >= triangle_buffer_size)
							{
								REPLAY_COUNT = 0;
							}
						REPLAY_COUNT += 1;
					}
				break;
			case 'S':
				if (DIRTY_MESH)
					{
						// works only in replay mode
						fprintf(stderr,"back stepping only work during second play (replay)\n");
					}
				else
					{
						NEXT_STEP = triangle_buffer_size / EXACTLY_N_STEPS;
						if (NEXT_STEP == 0) NEXT_STEP = 1;
						REPLAY_COUNT -= NEXT_STEP;
						if (REPLAY_COUNT < 0)
							{
								REPLAY_COUNT = 0;
							}
					}
				break;
			case 'P':
				DIRTY_MESH = 1;
			case 'p':
				if (DIRTY_MESH) {
					AnimationOn = !AnimationOn;
				}
				else {
					if (REPLAY_IT == 0) {
						if (REPLAY_COUNT >= triangle_buffer_size) {
							REPLAY_COUNT = 0;
						}
						NEXT_STEP = triangle_buffer_size / EXACTLY_N_STEPS;
						if (NEXT_STEP == 0) NEXT_STEP = 1;
						REPLAY_IT = 1;
					}
					else {
						REPLAY_IT = 0;
					}
				}
			case 's':
				if (DIRTY_MESH)
					{
						if (WorkingOn == 0)
							{
								vizBegin();
								WorkingOn = vizContinue();
							}
						else
							{
								WorkingOn = vizContinue();
							}
						if (WorkingOn == 0)
							{
								vizEnd();
								AnimationOn = 0;
								PrintToFileOn = 0;
							}
					}
				else
					{
						if (REPLAY_COUNT >= triangle_buffer_size)
							{
								REPLAY_COUNT = 0;
							}
						NEXT_STEP = triangle_buffer_size / EXACTLY_N_STEPS;
						if (NEXT_STEP == 0) NEXT_STEP = 1;
						REPLAY_COUNT += NEXT_STEP;
					}
				break;
			case 'L':
			case 'l':
				fprintf(stderr,"SHOW_LINES %d\n",!SHOW_LINES);
				if (SHOW_LINES)
					{
						SHOW_LINES = 0;
					}
				else
					{
						SHOW_LINES = 1;
						if (DIRTY_LINES)
							{
								fprintf(stderr,"computing lines ...\n");
								compute_lines();
								DIRTY_LINES = 0;
								DIRTY_VERTEX_SPAN = 1;
								SHOW_VERTEX_SPAN = 0;
							}
					}
				break;



			case '\'':
				if (SaveViewInfo("view.txt"))
					cerr<<"view information saved"<<endl;
				break;

			case ',':
				if (LoadViewInfo("view.txt"))
					cerr<<"view information loaded"<<endl;
				break;

			case '.':
				// have to do this in a separate thread!
				new thlib::Thread(rm_hack, this, 0);
				break;
			}

		Redraw();
	}


	static void* rm_hack(void *_viewer) {
		SMViewer *viewer = (SMViewer*)_viewer;

		for (int i=1; i<=5; i++) {
			char tmp[100];
			sprintf(tmp, "view.%d.txt", i);
			viewer->LoadViewInfo(tmp);
			viewer->doFull= true;
			viewer->Redraw();
			while (viewer->doFull) { }
			sprintf(tmp, "fullres.%d.ppm", i);
			rename("fullres.ppm", tmp);
		}

		return 0;
	}


	void render_triangle(GridVertex* vertex0, GridVertex* vertex1, GridVertex* vertex2)
	{
		float n[3];
		float v0[3];
		float v1[3];
		float v2[3];

		VecScalarDiv3fv(v0, vertex0->v, (float)vertex0->number);
		VecScalarDiv3fv(v1, vertex1->v, (float)vertex1->number);
		VecScalarDiv3fv(v2, vertex2->v, (float)vertex2->number);
		VecCcwNormal3fv(n, v0, v1, v2);

		glNormal3fv(n);
		glVertex3fv(v0);
		glVertex3fv(v1);
		glVertex3fv(v2);
	}

	void render_line(GridVertex* vertex0, GridVertex* vertex1)
	{
		float v0[3];
		float v1[3];

		VecScalarDiv3fv(v0, vertex0->v, (float)vertex0->number);
		VecScalarDiv3fv(v1, vertex1->v, (float)vertex1->number);

		glVertex3fv(v0);
		glVertex3fv(v1);
	}


	void BoundingSphere(int win, float &x, float &y, float &z, float &r) {
		if (smreader) {
			x = (smreader->bb_min_f[0] + smreader->bb_max_f[0]) / 2;
			y = (smreader->bb_min_f[1] + smreader->bb_max_f[1]) / 2;
			z = (smreader->bb_min_f[2] + smreader->bb_max_f[2]) / 2;

			float xl = (smreader->bb_min_f[0] - smreader->bb_max_f[0]) / 2;
			float yl = (smreader->bb_min_f[1] - smreader->bb_max_f[1]) / 2;
			float zl = (smreader->bb_min_f[2] - smreader->bb_max_f[2]) / 2;
			
			r = sqrt(xl*xl + yl*yl + zl*zl) * 1.2;
		} else {
			x = y = z = 0;
			r = 1;
		}
	}


	void Draw(int win)
	{

		if (doFull) {
			full_resolution_rendering();
			doFull = false;
			return;
		}

		SetUnsharpWidth(0);
		SetUnsharpLambda(4);

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
		// glViewport(0,0,WindowW,WindowH);

		// glMatrixMode(GL_PROJECTION);
		// glLoadIdentity();
		// gluPerspective(30.0f,(float)WindowW/WindowH,0.0625f,5.0f);

		glMatrixMode(GL_MODELVIEW);
		// glLoadIdentity();
		// gluLookAt(DistX,DistY,DistZ, DistX,DistY,0, 0,1,0);

		// glRotatef(Elevation,1,0,0);
		// glRotatef(Azimuth,0,1,0);

		// glTranslatef(boundingBoxTranslateX,boundingBoxTranslateY,boundingBoxTranslateZ*EXTRA_Z_SCALE);
		// glScalef(boundingBoxScale,boundingBoxScale,boundingBoxScale*EXTRA_Z_SCALE);

		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_NORMALIZE);

		if (RENDER_MODE != 3) {
			if (RENDER_MODE == 0) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
			else if (RENDER_MODE == 1) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			}
			else if (RENDER_MODE == 2) {
				glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
			}

			if (switch_colors) {
				glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_white);
				glMaterialfv(GL_BACK, GL_DIFFUSE, colours_light_blue);
			}
			else {
				glMaterialfv(GL_FRONT, GL_DIFFUSE, colours_light_blue);
				glMaterialfv(GL_BACK, GL_DIFFUSE, colours_white);
			}
  
			if (SHOW_LINES) {
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0f, 1.0f);
				glLineWidth(2.0f);
			}

			STRIPE_WHITE = 0;

			int rendered_triangles;

			if (DIRTY_MESH) {
				rendered_triangles = triangle_buffer_size;
			}
			else {
				if (REPLAY_COUNT > triangle_buffer_size) {
					rendered_triangles = triangle_buffer_size;
					REPLAY_IT = 0;
				}
				else {
					rendered_triangles = REPLAY_COUNT;
				}
			}

			// draw triangles

			if (triangle_buffer_size) {
				glBegin(GL_TRIANGLES);
				for (int i = 0; i < rendered_triangles; i++) {
					if (SHOW_VERTEX_SPAN) {
						int v_span_max = grid_vertex_buffer[triangle_buffer[3*i+0]].v_span;
						if (v_span_max < grid_vertex_buffer[triangle_buffer[3*i+1]].v_span) v_span_max = grid_vertex_buffer[triangle_buffer[3*i+1]].v_span;
						if (v_span_max < grid_vertex_buffer[triangle_buffer[3*i+2]].v_span) v_span_max = grid_vertex_buffer[triangle_buffer[3*i+2]].v_span;
						float red = 0.7f-(0.7f/65535*RED_SCALE*v_span_max);
						colours_red_white[1] = red;
						colours_red_white[2] = red;

						if (switch_colors) {
							glMaterialfv(GL_FRONT, GL_DIFFUSE,colours_red_white);
						}
						else {
							glMaterialfv(GL_BACK, GL_DIFFUSE, colours_red_white);
						}
					}
					else if (STREAM_COLORING) {
						if (STREAM_COLORING == 1) {
							if ((i % (triangle_buffer_size/NUMBER_STRIPES)) == 0) {
								if (STRIPE_WHITE) {
									glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,colours_stripe_white);
								}
								else {
									glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_stripe_black);
								}
								STRIPE_WHITE = !STRIPE_WHITE;
							}
						}
						else if (STREAM_COLORING == 2) {
							colours_slow_change[0] = 0.1f+0.8f*i/triangle_buffer_size;
							colours_slow_change[1] = colours_slow_change[0];
							colours_slow_change[2] = colours_slow_change[0];
							glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_slow_change);
						}
						else if (STREAM_COLORING == 3) {
							if (i < triangle_buffer_size/3) {
								colours_slow_change[0] = 0.1f+0.8f*i/(triangle_buffer_size/3);
								colours_slow_change[1] = 0.1f;
								colours_slow_change[2] = 0.1f;
							}
							else if (i < 2*(triangle_buffer_size/3)) {
								colours_slow_change[0] = 0.9f;
								colours_slow_change[1] = 0.1f+0.8f*(i-(triangle_buffer_size/3))/(triangle_buffer_size/3);
								colours_slow_change[2] = 0.1f;
							}
							else {
								colours_slow_change[0] = 0.9f;
								colours_slow_change[1] = 0.9f;
								colours_slow_change[2] = 0.1f+0.8f*(i-2*(triangle_buffer_size/3))/(triangle_buffer_size/3);
							}
							glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, colours_slow_change);
						}
					}
					render_triangle(&(grid_vertex_buffer[triangle_buffer[3*i+0]]), &(grid_vertex_buffer[triangle_buffer[3*i+1]]), &(grid_vertex_buffer[triangle_buffer[3*i+2]]));
				}
				glEnd();
			}
		}


		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_NORMALIZE);

		if (SHOW_LINES) {
			glColor3fv(colours_diffuse[0]);
			glBegin(GL_LINES);
			for (int i = 0; i < triangle_buffer_size; i++) {
				for (int j = 0; j < 3; j++) {
					if ((grid_vertex_buffer[triangle_buffer[3*i+j]].v_span == 3) && (grid_vertex_buffer[triangle_buffer[3*i+((j+1)%3)]].v_span == 3)) {
						render_line(&(grid_vertex_buffer[triangle_buffer[3*i+j]]), &(grid_vertex_buffer[triangle_buffer[3*i+((j+1)%3)]]));
					}
				}
			}
			glEnd();
			glDisable(GL_POLYGON_OFFSET_FILL);
		}


		glPointSize(2);
		glBegin(GL_POINTS);
		glColor3f(1,0,0);
		for (int v=0; v<boundary_verts.size(); v+=4) {
			if (boundary_verts[v+3] > 0) {
				glColor3f(1,0,0);
			} else {
				glColor3f(0,0,1);
			}
			glVertex3fv(&boundary_verts[v]);
		}
		glEnd();
		glPointSize(1);

		glDisable(GL_DEPTH_TEST);
	}
};


SMViewer *ggui;
void MyMenuFunc(int i) {
	ggui->MyMenuFunc(i);
}
 
int main(int argc, char *argv[])
{
	SMViewer gui(argc, argv);
	ggui = &gui;

	gui.Go();

	return 0;
}
