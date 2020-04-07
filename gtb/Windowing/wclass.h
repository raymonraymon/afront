
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


#ifndef __WCLASS_H_777
#define __WCLASS_H_777

#ifdef WIN32

#include <gtb/common.hpp>

/*
 * Q&D hack to open a window that I can draw openGL on
 * Used for offline rendering with glut applications.
 */
GTB_BEGIN_NAMESPACE

class GLWClass
{
public:
    GLWClass(const char* name, int x, int y, int width, int height);
    ~GLWClass();

    void SwapBuffers();
    void PrepareForRendering();
    void Show();
    void Hide();

protected:
    static int lastid;
    char classname[1000];

    HWND hWin;                  // Window handler
    HDC hDC;                    // Device context
    HGLRC hRC;                  // Rendering context

    static LRESULT WINAPI wproc( HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam );

    void EnableOpenGL(HWND hWnd, HDC * hDC, HGLRC * hRC);
    void DisableOpenGL(HWND hWnd, HDC hDC, HGLRC hRC);
};

GTB_END_NAMESPACE

#endif // WIN32

#endif // __WCLASS_H_777
