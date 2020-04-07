
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

#include "wclass.h"

GTB_BEGIN_NAMESPACE

int GLWClass::lastid = 1;

GLWClass::GLWClass(const char* name, int x, int y, int width, int height)
{
    sprintf(classname, "%s%d", name, lastid);
    ++lastid;

    WNDCLASS  wc;
    HINSTANCE hInstance = GetModuleHandle(NULL);
    /* Clear (important!) and then fill in the window class structure. */
    memset(&wc, 0, sizeof(WNDCLASS));
    wc.style         = CS_OWNDC;
    wc.lpfnWndProc   = (WNDPROC)wproc;
    wc.hInstance     = hInstance;
    wc.hIcon         = NULL; //LoadIcon(hInstance, "GLUT_ICON");
    wc.hCursor       = NULL; // LoadCursor(hInstance, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName  = NULL;
    wc.lpszClassName = classname;
    if(!RegisterClass(&wc)) {
        printf("Failed to register class\n");
    }

// Create the window
    hWin = CreateWindow(
        classname, name,
		WS_CAPTION | WS_POPUPWINDOW | WS_VISIBLE,
        x, y, width, height,
        NULL,
        NULL, 
        NULL,
        NULL);

    EnableOpenGL(hWin, &hDC, &hRC);
}


GLWClass::~GLWClass()
{
    DisableOpenGL(hWin, hDC, hRC);
    DestroyWindow(hWin);
    UnregisterClass(classname, GetModuleHandle(NULL));
}

void GLWClass::SwapBuffers()
{
    ::SwapBuffers( hDC );
}

void GLWClass::PrepareForRendering()
{
    if (wglGetCurrentDC() != hDC)
    {
        wglMakeCurrent(hDC, hRC);
    }
}


void GLWClass::EnableOpenGL(HWND hWnd, HDC * hDC, HGLRC * hRC)
{
	PIXELFORMATDESCRIPTOR pfd;
	int format;
	
	// get the device context (DC)
	*hDC = GetDC( hWnd );
	
	// set the pixel format for the DC
	ZeroMemory( &pfd, sizeof( pfd ) );
	pfd.nSize = sizeof( pfd );
	pfd.nVersion = 1;
	pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 24;
	pfd.cDepthBits = 16;
	pfd.iLayerType = PFD_MAIN_PLANE;
	format = ChoosePixelFormat( *hDC, &pfd );
	SetPixelFormat( *hDC, format, &pfd );
	
	// create and enable the render context (RC)
	*hRC = wglCreateContext( *hDC );
	wglMakeCurrent( *hDC, *hRC );
	
}

// Disable OpenGL

void GLWClass::DisableOpenGL(HWND hWnd, HDC hDC, HGLRC hRC)
{
//	wglMakeCurrent( NULL, NULL );
	wglDeleteContext( hRC );
	ReleaseDC( hWnd, hDC );
}

LRESULT WINAPI GLWClass::wproc( HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam )
{
	switch (msg)
	{
		
	case WM_CREATE:
		return 0;
		
	case WM_CLOSE:
		PostQuitMessage( 0 );
		return 0;
		
	case WM_DESTROY:
		return 0;
		
	case WM_KEYDOWN:
		switch ( wparam )
		{
			
		case VK_ESCAPE:
			PostQuitMessage(0);
			return 0;
			
		}
		return 0;
	
	default:
		return DefWindowProc( hwnd, msg, wparam, lparam );
	}
}

void GLWClass::Show()
{
    ShowWindow(hWin, SW_SHOW);
}

void GLWClass::Hide()
{
    ShowWindow(hWin, SW_HIDE);
}

GTB_END_NAMESPACE
