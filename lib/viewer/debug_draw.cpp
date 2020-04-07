
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
#include "debug_draw.h"
#include <viewer/viewer.h>

using namespace std;
using namespace gtb;

// polygons
DbgPoly DbgPoly::wins[2];
void DbgPoly::madd(const std::vector<gtb::tPoint3<float> > &p, const std::vector<gtb::tPoint3<float> > &c){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.push_back(p); colors.push_back(c);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPoly::translate(const gtb::tVector3<float> &t) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    for (unsigned l=0; l<points.size(); l++) {
	for (unsigned p=0; p<points[l].size(); p++) {
	    points[l][p] = points[l][p] + t;
	}
    }
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPoly::clear(){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.resize(0); colors.resize(0);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPoly::draw() const {
    for (unsigned i=0; i<points.size(); i++) {
	(points[i][1]-points[i][0]).cross(points[i][2]-points[i][0]).normalized().load_as_normal();
	glBegin(GL_POLYGON);
	for (unsigned p=0; p<points[i].size(); p++) {
	    glColor3f(colors[i][p][0], colors[i][p][1], colors[i][p][2]);
	    points[i][p].load();
	}
	glEnd();
    }
}



// poly lines
DbgPLines DbgPLines::wins[2];
void DbgPLines::madd(const std::vector<gtb::tPoint3<float> > &a, int c, int w){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    lines.push_back(a);		color.push_back(c);		width.push_back(w);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPLines::translate(const gtb::tVector3<float> &t) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    for (unsigned l=0; l<lines.size(); l++) {
	for (unsigned p=0; p<lines[l].size(); p++) {
	    lines[l][p] = lines[l][p] + t;
	}
    }
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPLines::clear(){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    lines.resize(0); color.resize(0); width.resize(0);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPLines::draw() const {
    for (unsigned i=0; i<lines.size(); i++) {
	switch (color[i]) {
	case 0:	glColor3f(1,0,0);	break;
	case 1:	glColor3f(0,1,0);	break;
	case 2:	glColor3f(0,0,1);	break;
	case 3:	glColor3f(0,0,0);	break;
	}

	glLineWidth(width[i]);
	glBegin(GL_LINE_STRIP);
	for (unsigned p=0; p<lines[i].size(); p++) {
	    lines[i][p].load();
	}
	glEnd();
    }
	glLineWidth(1);
}


// points
DbgPoints DbgPoints::wins[2];
void DbgPoints::madd(const gtb::tPoint3<float> &p, float r, float g, float b) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.push_back(p);
    color.push_back(r); color.push_back(g); color.push_back(b);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPoints::translate(const gtb::tVector3<float> &t) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    for (unsigned p=0; p<points.size(); p++) {
	points[p] = points[p] + t;
    }
    Viewer::GetViewer()->GetCriticalSection().leave();
}

void DbgPoints::clear(){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.resize(0);
    color.resize(0);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgPoints::draw() const {
    if (!Viewer::GetViewer()) return;
    glPointSize(4);
    glBegin(GL_POINTS);

    for (unsigned i=0; i<points.size(); i++) {
	glColor3f(color[i*3+0], color[i*3+1], color[i*3+2]);
	points[i].load();
    }

    glEnd();
    glPointSize(1);
}


// oriented points
DbgOPoints DbgOPoints::wins[2];
void DbgOPoints::madd(const gtb::tPoint3<float> &p, const gtb::tVector3<float> &n, float r, float g, float b) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.push_back(p);
    normals.push_back(n);
    color.push_back(r); color.push_back(g); color.push_back(b);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgOPoints::translate(const gtb::tVector3<float> &t) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    for (unsigned p=0; p<points.size(); p++) {
	points[p] = points[p] + t;
    }
    Viewer::GetViewer()->GetCriticalSection().leave();
}

void DbgOPoints::clear(){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.resize(0);
    normals.resize(0);
    color.resize(0);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgOPoints::draw() const {
    if (!Viewer::GetViewer()) return;
    glPointSize(4);
    glBegin(GL_POINTS);

    for (unsigned i=0; i<points.size(); i++) {
	glColor3f(color[i*3+0], color[i*3+1], color[i*3+2]);
	normals[i].load_as_normal();
	points[i].load();
    }

    glEnd();
    glPointSize(1);
}


// spheres
DbgSpheres DbgSpheres::wins[2];
void DbgSpheres::madd(const gtb::tPoint3<float> &p, float rad, float r, float g, float b, float a) {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.push_back(p);
    radius.push_back(rad);
    color.push_back(r); color.push_back(g); color.push_back(b); color.push_back(a);
    Viewer::GetViewer()->GetCriticalSection().leave();
}

void DbgSpheres::clear(){
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    points.resize(0);
    radius.resize(0);
    color.resize(0);
    Viewer::GetViewer()->GetCriticalSection().leave();
}
void DbgSpheres::draw() const {
    static GLUquadricObj* gluobj = gluNewQuadric();

    glPushAttrib(GL_ENABLE_BIT);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (unsigned i=0; i<points.size(); i++) {
	glColor4f(color[i*4+0], color[i*4+1], color[i*4+2], color[i*4+3]);

	glPushMatrix();
	glTranslated(points[i][0], points[i][1], points[i][2]);
	glScalef(radius[i], radius[i], radius[i]);
	gluSphere(gluobj, 1, 10, 10);
	glPopMatrix();
    }

    glPopAttrib();
}




void dbgClear() {
    if (!Viewer::GetViewer()) return;
    Viewer::GetViewer()->GetCriticalSection().enter();
    DbgPLines::clear(0);
    DbgPLines::clear(1);
    DbgPoints::clear(0);
    DbgPoints::clear(1);
    DbgOPoints::clear(0);
    DbgOPoints::clear(1);
    DbgSpheres::clear(0);
    DbgSpheres::clear(1);
    DbgPoly::clear(0);
    DbgPoly::clear(1);
    Viewer::GetViewer()->GetCriticalSection().leave();
}


