
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



#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>

#include <gtb/graphics/ogltools.h>

#include <math.h>

#include "viewer.h"
#include "trackball.h"


void *Viewer::start_stub(void *arg)
{
    theViewer->Go();
    return 0;
}


void Viewer::Start() {
	_thread = new thlib::Thread(start_stub, NULL, 0);
	while (!glut_has_initialized); // Busywaiting, UGH
}




void Viewer::Go() {
    _critical_section.enter();

	int argc=1;
	char *argv[1];
	char blah[100] = "blah";
	argv[0]=blah;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);

	for (unsigned i=0; i<windows.size(); i++) {
	    MakeWindow(i);
	}
	glut_has_initialized = true;
	_critical_section.leave();
	

    glutMainLoop();	// never returns!
}

void *Viewer::cam_stub(void *arg)
{
	theViewer->CamSyncLoop();
	return 0;
}

void Viewer::CamSyncLoop() {

	while (1) {

		float pan[3];
		float d;
		float quat[4];
		bool o;

		if (ReadViewInfo(stdin, pan, d, quat, o)) {
			for (unsigned i=0; i<windows.size(); i++) {
				if (windows[i]->sync_cam_in)
					LoadViewInfo(pan, d, quat, o, i);
			}
		}

		thlib::Thread::yield();
	}
}



Viewer* Viewer::theViewer = NULL;
volatile bool Viewer::glut_has_initialized = false;

Viewer::Viewer() :
_critical_section(),
_thread(NULL),
cam_thread(NULL),
input_cs(),
last_key(-1),
last_key_win(-1),
last_click_win(-1),
last_click_button(-1),
key_send(0),
key_send_window(-1)
{
	if (theViewer) {
		std::cerr<<"can only create one viewer, exiting!"<<std::endl;
		exit(-1);
	}
	theViewer=this;
}

Viewer::~Viewer() {
#ifndef WIN32
    if (_thread != NULL) {
	pthread_kill(*(pthread_t*)_thread, SIGKILL);

    }
#endif
}


int Viewer::AddWindow(const char *title, int width, int height, int xpos, int ypos) {
	windows.push_back(new Window(title, xpos, ypos, width, height));
	return ((int)windows.size()-1);
}

void Viewer::MakeWindow(int win)
{
    Window &w = *windows[win];
    
    glutInitWindowSize(w.iwidth, w.iheight);
    glutInitWindowPosition(w.ixpos, w.iypos);
    w.glut_id = glutCreateWindow(w.title);

    // register callbacks
    glutDisplayFunc(s_gc_display);
    glutReshapeFunc(s_gc_reshape);
    glutKeyboardFunc(s_gc_keyboard);
    glutMouseFunc(s_gc_mouse);
    glutMotionFunc(s_gc_motion);
    glutPassiveMotionFunc(s_gc_passive_motion);
    glutIdleFunc(s_gc_idle);

	glutSetCursor(GLUT_CURSOR_CROSSHAIR);
    Init(win);
}


void Viewer::Redraw(int win) {
	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}


void Viewer::ToggleImmediateMode(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->immediatemode = !windows[i]->immediatemode;
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}


void Viewer::IncUnsharpWidth(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_width++;
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}


void Viewer::DecUnsharpWidth(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_width = std::max(1, windows[i]->unsharp_width-1);
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}

void Viewer::SetUnsharpWidth(int w, int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_width = w;
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}

void Viewer::IncUnsharpLambda(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_lambda++;
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}


void Viewer::DecUnsharpLambda(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_lambda = std::max(1, windows[i]->unsharp_lambda-1);
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}

void Viewer::SetUnsharpLambda(int w, int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->unsharp_lambda = w;
		windows[i]->redraw=true;
		windows[i]->cs.leave();
	}
}


void Viewer::ApplyUnsharp(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int _w=begin; _w<end; _w++) {
		Window &w = *windows[_w];
		w.cs.enter();

		// get the color and depth buffers
		float *color = new float[w.width*w.height*4];

		float *depth = new float[w.width*w.height];
		
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, w.width, w.height, GL_DEPTH_COMPONENT, GL_FLOAT, (void*)depth);
		glReadPixels(0, 0, w.width, w.height, GL_RGBA, GL_FLOAT, (void*)color);

		float n = fabs((-w.pan_z+w.dolly))/1000;
		float f = (-w.pan_z + w.dolly+2)*2;

		// unproject the z's
		for (int y=0; y<w.height; y++) {
			for (int x=0; x<w.width; x++) {
				GLdouble objx, objy, objz;
				if (depth[y*w.width+x] == 1) {
					depth[y*w.width+x] = 0;
				} else {
					float z = (2 * n) / (f + n - depth[y*w.width+x] * (f - n));
					depth[y*w.width+x] = - (f*n) / ((depth[y*w.width+x]/* / s*/) * (f-n) - f);
				}
			}
		}

		// rescale the depths to be in [0,1]
		float min=1000000;
		float max=-1000000;
		for (int i=0; i<w.width*w.height; i++) {
			if (depth[i] != 0) {
				min = std::min(min, depth[i]);
				max = std::max(max, depth[i]);
			}
		}


		for (int i=0; i<w.width*w.height; i++) {
			if (depth[i] != 0) {
				depth[i] = (depth[i]-min)/(max-min);
			}
		}



		// blur the depth buffer
		float *depth_b = new float[w.width*w.height];
		float *depth_t = new float[w.width*w.height]; // transposed!

		// first in x
		for (int y=0; y<w.height; y++) {

			for (int x=0; x<w.width; x++) {
				depth_t[x*w.height+y] = 0;
				float weight = 0;

				for (int xw=std::max(0,x-w.unsharp_width);
					 xw<=std::min(w.width-1,x+w.unsharp_width);
					 xw++) {

					float thisw = 1;
					float thisd = depth[y*w.width + xw];

					if (thisd != 0) {
						depth_t[x*w.height+y] += thisw * thisd;
						weight += thisw;
					}
				}

				if (weight!=0 && depth[y*w.width+x]!=0)
					depth_t[x*w.height+y] /= weight;
				else
					depth_t[x*w.height+y] = depth[y*w.width+x];
			}
		}


		// now in y, but from the transposed
		for (int x=0; x<w.width; x++) {

			for (int y=0; y<w.height; y++) {
				depth_b[y*w.width+x] = 0;
				float weight = 0;

				for (int yw=std::max(0,y-w.unsharp_width);
					 yw<=std::min(w.height-1,y+w.unsharp_width);
					 yw++) {

					float thisw = 1;
					float thisd = depth_t[x*w.height + yw];

					if (thisd != 0) {
						depth_b[y*w.width + x] += thisw * thisd;
						weight += thisw;
					}
				}

				if (weight!=0 && depth[y*w.width+x]!=0)
					depth_b[y*w.width+x] /= weight;
				else
					depth_b[y*w.width+x] = depth[y*w.width+x];
			}
		}
		delete [] depth_t;



		// compute the new colors
		int tw = 1;
		while (tw < w.width) tw<<=1;

		int th = 1;
		while (th < w.height) th<<=1;

		float *ncolor = new float[tw*th*4];
		for (int i=0; i<tw*th*4; i++) ncolor[i] = 0;

		for (int y=0; y<w.height; y++) {
			for (int x=0; x<w.width; x++) {
				float dd = depth_b[y*w.width+x] - depth[y*w.width+x];
				if (dd >= 0) {
					ncolor[(y*tw + x)*4 + 0] = color[(y*w.width + x)*4 + 0];
					ncolor[(y*tw + x)*4 + 1] = color[(y*w.width + x)*4 + 1];
					ncolor[(y*tw + x)*4 + 2] = color[(y*w.width + x)*4 + 2];
				} else {

					if (dd<0) dd=-dd;
					dd = exp(-pow(dd*w.unsharp_lambda,2))*0.25+0.75;

					ncolor[(y*tw + x)*4 + 0] = color[(y*w.width + x)*4 + 0] * (dd);
					ncolor[(y*tw + x)*4 + 1] = color[(y*w.width + x)*4 + 1] * (dd);
					ncolor[(y*tw + x)*4 + 2] = color[(y*w.width + x)*4 + 2] * (dd);
				}
			}
		}


		// render the new colors back
		GLuint texture;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		
		glTexImage2D(GL_TEXTURE_2D,
					 0,
					 GL_RGB,
					 tw,
					 th,
					 0,
					 GL_RGBA,
					 GL_FLOAT,
					 ncolor);


		glEnable(GL_TEXTURE_2D);
		glDisable(GL_DEPTH_TEST);
		glDepthMask(0);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0,w.width,0,w.height);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		glBegin(GL_QUADS);

		glColor3f(1,1,1);

		glTexCoord2f(0,0);
		glVertex2f(0,0);

		glTexCoord2f((float)w.width/tw, 0);
		glVertex2f(w.width,0);

		glTexCoord2f((float)w.width/tw, (float)w.height/th);
		glVertex2f(w.width,w.height);

		glTexCoord2f(0,(float)w.height/th);
		glVertex2f(0,w.height);

		glEnd();
		

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_TEXTURE_2D);
		glDepthMask(1);

		glDeleteTextures(1,&texture);


		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();


		delete [] depth;
		delete [] depth_b;

		delete [] ncolor;
		delete [] color;



		w.cs.leave();
	}
}


void Viewer::ToggleOrtho(int win) {

	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->ortho = !windows[i]->ortho;
		windows[i]->redraw=true;
		windows[i]->cs.leave();

		if (windows[i]->sync_cam_out)
			SaveViewInfo(stdout, i);
	}
}


void Viewer::WaitForKey(int key, int win) {

	input_cs.enter();
	last_key=-1;
	input_cs.leave();


	while (1) {

		thlib::Thread::yield();
		input_cs.enter();
		if (last_key >= 0) {	// key was hit
			if ((win<0 || win==last_key_win) && // want any win, or the win matches
				(key<=0 || key==last_key)) {		// want any key, or the key matches
				input_cs.leave();
				Redraw(win);
			    break;
			}
		}
		input_cs.leave();
	}

}


void Viewer::WaitForKeyOrClick(int key, int kwin, int button, int bwin) {

	input_cs.enter();
	last_key=-1;
	input_cs.leave();


	while (1) {

		thlib::Thread::yield();
		input_cs.enter();
		if (last_key >= 0) {	// key was hit
			if ((kwin<0 || kwin==last_key_win) && // want any win, or the win matches
				(key<=0 || key==last_key)) {		// want any key, or the key matches
				input_cs.leave();
				Redraw(kwin);
			    break;
			}
		}

		if (last_click_button >= 0) { // mouse was clicked
			if ((bwin<0 || bwin==last_click_win) &&
				(button<0 || button==last_click_button)) {
				input_cs.leave();
				Redraw(bwin);
			    break;
			}
		}


		input_cs.leave();
	}

}


int Viewer::LastKey() {

	input_cs.enter();

	int ret = last_key;
	last_key = -1;

	input_cs.leave();

	return ret;
}


void Viewer::SetCamSyncIn(int win) {
	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->sync_cam_in=true;
		windows[i]->cs.leave();
	}

	if (!cam_thread) {
		cam_thread = new thlib::Thread(cam_stub, NULL, 0);
	}

}

void Viewer::SetCamSyncOut(int win) {
	int begin=(win<0) ? 0 : win;
	int end = (win<0) ? (int)windows.size() : win+1;
	
	for (int i=begin; i<end; i++) {
		windows[i]->cs.enter();
		windows[i]->sync_cam_out=true;
		windows[i]->cs.leave();
	}
}


bool Viewer::SaveViewInfo(const char *filename, int win) {
	FILE *f = fopen(filename, "w");
	if (!f) return false;
	SaveViewInfo(f, win);
	fclose(f);
	return true;
}

void Viewer::SaveViewInfo(FILE *f, int win) {

	Window &w = *windows[win];
	w.cs.enter();
	fprintf(f, "pan: %f %f %f\n", w.pan_x, w.pan_y, w.pan_z);
	fprintf(f, "dolly: %f\n", w.dolly);
	fprintf(f, "quat: %f %f %f %f\n", w.quat[0], w.quat[1], w.quat[2], w.quat[3]);
	fprintf(f, "ortho: %i\n", w.ortho);
	fflush(f);
	w.cs.leave();
}


bool Viewer::LoadViewInfo(const char *filename, int win) {

	FILE *f = fopen(filename, "r");
	if (!f) return false;

	if (!LoadViewInfo(f, win)) {
		fclose(f);
		return false;
	}
	fclose(f);
	return true;
}


bool Viewer::LoadViewInfo(FILE *f, int win) {

	float pan[3];
	float d;
	float quat[4];
	bool o;

	if (!ReadViewInfo(f, pan, d, quat, o)) return false;
	LoadViewInfo(pan, d, quat, o, win);

	return true;
}

bool Viewer::ReadViewInfo(FILE *f, float pan[3], float &dolly, float quat[4], bool &ortho) {
	char line[1024];

	if (!fgets(line, 1000, f)) return false;
	if (3 != sscanf(line, "pan: %f %f %f", &pan[0], &pan[1], &pan[2])) {
		std::cerr<<"ReadViewInfo: couldn't parse camera sting: "<<line<<std::endl;
		return false;
	}

	if (!fgets(line, 1000, f)) return false;
	if (1 != sscanf(line, "dolly: %f", &dolly)) {
		std::cerr<<"ReadViewInfo: couldn't parse camera sting: "<<line<<std::endl;
		return false;
	}

	if (!fgets(line, 1000, f)) return false;
	if (4 != sscanf(line, "quat: %f %f %f %f", &quat[0], &quat[1], &quat[2], &quat[3])) {
		std::cerr<<"ReadViewInfo: couldn't parse camera sting: "<<line<<std::endl;
		return false;
	}

	if (!fgets(line, 1000, f)) return false;
	if (1 != sscanf(line, "ortho: %i", &ortho)) {
		std::cerr<<"ReadViewInfo: couldn't parse camera sting: "<<line<<std::endl;
		return false;
	}

	return true;
}


void Viewer::LoadViewInfo(float pan[3], float dolly, float quat[4], bool ortho, int win) {
	Window &w = *windows[win];
	w.cs.enter();

	w.pan_x = pan[0];
	w.pan_y = pan[1];
	w.pan_z = pan[2];

	w.dolly = dolly;

	w.quat[0] = quat[0];
	w.quat[1] = quat[1];
	w.quat[2] = quat[2];
	w.quat[3] = quat[3];

	w.ortho = ortho;

	w.redraw = true;
	w.cs.leave();
}


void Viewer::ViewAll(int win) {
    int begin = (win<0) ? 0 : win;
    int end = (win<0) ? (int)windows.size() : win+1;

	for (int i=begin; i<end; ++i)  {
		Window &w = *windows[i];

		float x,y,z,r;
		BoundingSphere(i, x,y,z,r);

		w.cs.enter();
		w.dolly = 2*r;
		w.pan_x = -x;
		w.pan_y = -y;
		w.pan_z = -z;
		w.quat[3]=1;
		w.quat[0]=w.quat[1]=w.quat[2] = 0;
		w.cs.leave();

		Redraw(i);

		if (w.sync_cam_out)
			SaveViewInfo(stdout, i);
	}

}

void Viewer::LookAt(int win, float x, float y, float z) {
	Window &w = *windows[win];
	w.cs.enter();
	w.pan_x = -x;
	w.pan_y = -y;
	w.pan_z = -z;
	w.cs.leave();

	if (w.sync_cam_out)
		SaveViewInfo(stdout, win);
}



void Viewer::ScreenShot(const char *filename, bool front, int win) {
	std::cerr<<"saving screenshot: "<<filename<<std::endl;

	Window &w = *windows[win];
	w.cs.enter();

	unsigned char *buff = new unsigned char[w.width*w.height*3];

	if (front) {
		glReadBuffer(GL_FRONT);
	} else {
		glReadBuffer(GL_BACK);
	}
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, w.width, w.height, GL_RGB, GL_UNSIGNED_BYTE, buff);

	w.cs.leave();


	FILE *f = fopen(filename, "w");
	if (!f) {
		std::cerr<<"couldn't open file: "<<filename<<std::endl;
	} else {
		fprintf(f, "P3 %d %d\n 255\n", w.width, w.height);

		for (int y=w.height-1; y>=0; y--) {
			for (int x=0; x<w.width; x++) {
				int i = (y*w.width + x)*3;
				fprintf(f, "%d %d %d\n", buff[i], buff[i+1], buff[i+2]);
			}
		}

		fclose(f);
	}

	delete [] buff;
}


void Viewer::SetMessage(int win, const char *msg) {
	Window &w = *windows[win];

	w.cs.enter();
	if (!msg) {
		w.status_msg[0] = 0;
	} else {
		strcpy(w.status_msg, msg);
	}
	w.cs.leave();
}


int Viewer::GetCurrentWindowIndex() {
	int id = glutGetWindow();

	for (int i=0; i<(int)windows.size(); i++) {
		if (windows[i]->glut_id == id)
			return i;
	}

	std::cerr<<"unknown current window, exiting!!"<<std::endl;
	exit(-1);
	return -1;
}


void Viewer::s_gc_display() {
	theViewer->gc_display();
}
void Viewer::gc_display() {

	if (key_send != 0) {
		int begin=(key_send_window<0) ? 0 : key_send_window;
		int end = (key_send_window<0) ? (int)windows.size() : key_send_window+1;
	
		for (int i=begin; i<end; i++) {
			Key(i, key_send);
		}

		key_send = 0;
	}


	int win = GetCurrentWindowIndex();
	Window &w = *windows[win];

	w.cs.enter();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if (w.ortho) {
		glOrtho(-(w.dolly), (w.dolly), -(GLfloat)w.height / (GLfloat)w.width*(w.dolly), (GLfloat)w.height / (GLfloat)w.width*(w.dolly), fabs((-w.pan_z+w.dolly))/100, (-w.pan_z + w.dolly+2)*2);
	} else {
		gluPerspective(45, (GLfloat)w.width / (GLfloat)w.height, fabs((-w.pan_z+w.dolly))/1000, (-w.pan_z + w.dolly+2)*2);
	}


    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	float light_dir[4] = { 0.2f, 0.2f, 1, 0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_dir);

	glTranslatef(0,0,-w.dolly);
  	GLfloat m[4][4];
	build_rotmatrix(m, w.quat);
	glMultMatrixf(&m[0][0]);
	glTranslatef(w.pan_x, w.pan_y, w.pan_z);


// 	if (w.immediatemode) {
	Draw(win);
	w.redraw=false;
// 	} else {
// 		if (w.redraw || w.display_list<0) {
// 			if (w.display_list<0)
// 				w.display_list = glGenLists(1);
// 			w.redraw=false;

// 			glNewList(w.display_list, GL_COMPILE_AND_EXECUTE);
// 			Draw(win);
// 			glEndList();
// 		} else {
// 			glCallList(w.display_list);
// 		}
// 	}



	if (!w.dollying && !w.tracking && !w.panning && w.unsharp_width > 1) {
		ApplyUnsharp(win);
	}


	// draw the status message
	glColor3f(0,0,0);
	glDisable(GL_LIGHTING);

	{

		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0,w.width,0,w.height);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glRasterPos3f(35,75,0);

		for(char *c=w.status_msg; (*c && *c!='\n'); c++)
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*c);

		glEnable(GL_DEPTH_TEST);


		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}


	w.cs.leave();

	glutSwapBuffers();
}


void Viewer::s_gc_reshape(int width, int height) {
	theViewer->gc_reshape(width, height);
}
void Viewer::gc_reshape(int width, int height) {
	int win = GetCurrentWindowIndex();
	windows[win]->cs.enter();
	windows[win]->width = width;
	windows[win]->height= height;
	windows[win]->cs.leave();
	glViewport(0, 0, width, height);
	//	glutPostRedisplay();
}


void Viewer::s_gc_keyboard(unsigned char key, int x, int y) {
	theViewer->gc_keyboard(key, x, y);
}
void Viewer::gc_keyboard(unsigned char key, int x, int y) {
	int win = GetCurrentWindowIndex();
	Key(win, key);
	input_cs.enter();
	last_key = key;
	last_key_win = win;
	input_cs.leave();
	//	glutPostRedisplay();
}


void Viewer::s_gc_mouse(int button, int state, int x, int y) {
	theViewer->gc_mouse(button, state, x, y);
}
void Viewer::gc_mouse(int button, int state, int x, int y) {

	int win = GetCurrentWindowIndex();
	Window &w = *windows[win];

	w.cs.enter();

	w.mouse_x = x;
	w.mouse_y = y;

	double objx, objy, objz;
	float depth; 
	{
		GLdouble model[16], proj[16];
		GLint view[4];

		glGetDoublev(GL_MODELVIEW_MATRIX, model);
		glGetDoublev(GL_PROJECTION_MATRIX, proj);
		glGetIntegerv(GL_VIEWPORT, view);

		glReadBuffer(GL_FRONT);
		glReadPixels(x,w.height-y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,(void*)&depth);
		gluUnProject((GLdouble)x, (GLdouble)(w.height-y), depth, model, proj, view, &objx, &objy, &objz);
		if(depth>.9999) objz = 1e34;
	}


	// all we handle here is the left button
	if (button == GLUT_LEFT_BUTTON) {
		if (state == GLUT_DOWN) {
			w.pan_depth = depth;
//			std::cerr<<depth<<std::endl;
			if (w.pan_depth==1) w.pan_depth=0.999;
			if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
				w.dollying = true;
			else if (glutGetModifiers() == 0)
				w.tracking = true;
			else if (glutGetModifiers()==GLUT_ACTIVE_CTRL || glutGetModifiers()==GLUT_ACTIVE_ALT)
				w.panning = true;
		}

		if (state == GLUT_UP) {
			w.tracking=w.panning=w.dollying=false;
		}
	}


	if (state==GLUT_DOWN)
		w.mouse_button = button;
	else
		w.mouse_button = -1;

	w.cs.leave();

	if (state == GLUT_DOWN) {
		input_cs.enter();
		last_click_point[0] = objx;
		last_click_point[1] = objy;
		last_click_point[2] = objz;
		last_click_win = win;
		last_click_button = button;
		input_cs.leave();

	}

	// pass all clicks on
	Mouse(win, button, (state==GLUT_DOWN)?-1:1, x, y, depth, objx, objy, objz);

	if (w.sync_cam_out)
		SaveViewInfo(stdout, win);

	//	glutPostRedisplay();
}


void Viewer::s_gc_motion(int x, int y) {
	theViewer->gc_motion(x,y);
}
float Viewer::get_object_space(int x, int y, double &objx, double &objy, double &objz) 
{
	int win = GetCurrentWindowIndex();
	Window &w = *windows[win];
	w.cs.enter();
	float depth; 
	{
		GLdouble model[16], proj[16];
		GLint view[4];

		glGetDoublev(GL_MODELVIEW_MATRIX, model);
		glGetDoublev(GL_PROJECTION_MATRIX, proj);
		glGetIntegerv(GL_VIEWPORT, view);

		glReadBuffer(GL_FRONT);
		glReadPixels(x,w.height-y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,(void*)&depth);
		gluUnProject((GLdouble)x, (GLdouble)(w.height-y), depth, model, proj, view, &objx, &objy, &objz);
		if(depth>.9999) objz = 1e34;
	}
	w.cs.leave();
	return depth;
}

void Viewer::s_gc_passive_motion(int x, int y)
{
    theViewer->gc_passive_motion(x, y);
}

void Viewer::gc_passive_motion(int x, int y)
{
	int win = GetCurrentWindowIndex();

	double objx, objy, objz;
	get_object_space(x, y, objx, objy, objz);

	PassiveMotion(win, x, y, objx, objy, objz);
}


void Viewer::gc_motion(int x, int y) {

	int win = GetCurrentWindowIndex();
	Window &w = *windows[win];

	double objx, objy, objz;
	float depth = get_object_space(x, y, objx, objy, objz);

	if (w.dollying) {
		w.dolly *= exp((double)(w.mouse_y-y) * 2. / w.height);
	}


	w.cs.enter();

	if (w.panning) {
		GLdouble model[16], proj[16];
		GLint view[4];
		glGetDoublev(GL_MODELVIEW_MATRIX, model);
                glGetDoublev(GL_PROJECTION_MATRIX, proj);
                glGetIntegerv(GL_VIEWPORT, view);

		GLdouble ox, oy, oz;
		GLdouble nx, ny, nz;
		gluUnProject((GLdouble)w.mouse_x, (GLdouble)(w.height-w.mouse_y), w.pan_depth, model, proj, view, &ox, &oy, &oz);
		gluUnProject((GLdouble)x,		  (GLdouble)(w.height-y),		  w.pan_depth, model, proj, view, &nx, &ny, &nz);
		w.pan_x += nx-ox;		w.pan_y += ny-oy;		w.pan_z += nz-oz;
	}


	if (w.tracking) {

		float quat[4];
		trackball(quat, (2.0 * w.mouse_x - w.width) / w.width,
						(w.height - 2.0 * w.mouse_y) / w.height,
						(2.0 * x - w.width) / w.width,
						(w.height- 2.0 * y) / w.height);

		add_quats(quat, w.quat, w.quat);
	}


	w.mouse_x = x;
	w.mouse_y = y;



	// pass all it on
	Mouse(win, w.mouse_button, 0, x, y, depth, objx, objy, objz);

	w.cs.leave();

	if (w.sync_cam_out)
		SaveViewInfo(stdout, win);


	//	glutPostRedisplay();
}


void Viewer::s_gc_idle() {
	theViewer->gc_idle();
}
void Viewer::gc_idle() {

	bool didsomething=false;

	for (unsigned i=0; i<windows.size(); i++) {
		
		 Idle();

		 windows[i]->cs.enter();
		 if (windows[i]->redraw) {
			 didsomething=true;
			 glutPostWindowRedisplay(windows[i]->glut_id);
		 }
		 windows[i]->cs.leave();
	 }

	 if (!didsomething)
		 thlib::Thread::yield();
}

