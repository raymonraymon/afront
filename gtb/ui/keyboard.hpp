
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


#ifndef GTB_KEYBOARD_INCLUDED
#define GTB_KEYBOARD_INCLUDED

#include <gtb/common.hpp>


#ifdef WIN32

GTB_BEGIN_NAMESPACE

class keyboard {
public:
	enum t_modifiers {key_ctrl, key_shift, key_alt};

	// Hack to work around bug in GLUT, which does not generate
	// KeyPress or KeyRelease events fot the modifiers.
	static bool is_key_pressed(t_modifiers key);
	static bool is_shift_pressed();
	static bool is_control_pressed();
	static bool is_alt_pressed();
};

GTB_END_NAMESPACE

#else // WIN32

#include <X11/Xlib.h>
#include <X11/keysym.h>

GTB_BEGIN_NAMESPACE

class keyboard {
public:
	// Hack to work around bug in GLUT, which does not generate
	// KeyPress or KeyRelease events fot the modifiers.
	static bool is_key_pressed(int keysym);
	static bool is_shift_pressed();
	static bool is_control_pressed();
	static bool is_alt_pressed();
	static void print_keys_pressed();

protected:
	static Display *init_display();
	static Display *_disp;
};

GTB_END_NAMESPACE

#endif // WIN32


#endif // GTB_KEYBOARD_INCLUDED
