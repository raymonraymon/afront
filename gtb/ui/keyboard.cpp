
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
#ifndef WIN32
#include <gtb/ui/keyboard.hpp>
#include <assert.h>
#include <stdio.h>
//#include <stdlib.h>
#include <X11/keysym.h>
#endif // WIN32

GTB_BEGIN_NAMESPACE

#ifndef WIN32


Display *keyboard::_disp = keyboard::init_display();


Display *keyboard::init_display()
{
	Display *disp = XOpenDisplay(NULL);
	assert(disp != NULL);
	return disp;
}


bool keyboard::is_key_pressed(int keysym)
{
	static char keys_return[32];
	assert(_disp != NULL);
	XQueryKeymap(_disp, keys_return);
	int keycode = XKeysymToKeycode(_disp, keysym);
	int i = keycode / 8;
	int j = keycode % 8;
	return keys_return[i] & (1 << j);
}


bool keyboard::is_shift_pressed()
{
	return (is_key_pressed(XK_Shift_L) ||
		is_key_pressed(XK_Shift_R));
}


bool keyboard::is_control_pressed()
{
	return (is_key_pressed(XK_Control_L) ||
		is_key_pressed(XK_Control_R));
}


bool keyboard::is_alt_pressed()
{
	return (is_key_pressed(XK_Alt_L) ||
		is_key_pressed(XK_Alt_R) ||
		is_key_pressed(XK_Mode_switch)); // HACK
}


void keyboard::print_keys_pressed()
{
	static char keys_return[32];
	assert(_disp != NULL);
	XQueryKeymap(_disp, keys_return);
	bool any_pressed = false;
	for (unsigned i = 0; i < 32; i++) {
		for (unsigned j = 0; j < 8; j++) {
			if (keys_return[i] & (1 << j)) {
				any_pressed = true;
				int keycode = i * 8 + j;
				KeySym keysym =
				    XKeycodeToKeysym(_disp, keycode, 0);
				fprintf(stderr, "%s ",
					XKeysymToString(keysym));
			}
		}
	}
	if (any_pressed) {
		fprintf(stderr, "\n");
	}
}


#else // WIN32


bool keyboard::is_key_pressed(t_modifiers key)
{
	switch (key) {
	case key_ctrl:
		return GetKeyState(VK_LCONTROL) < 0;
	case key_shift:
		return GetKeyState(VK_LSHIFT) < 0;
	case key_alt:
		return GetKeyState(VK_LMENU) < 0;
	default:
		fprintf(stderr, "unrecognized key\n");
		exit(EXIT_FAILURE);
        return 0;
	}
}


bool keyboard::is_shift_pressed()
{
	return is_key_pressed(key_shift);
}


bool keyboard::is_control_pressed()
{
	return is_key_pressed(key_ctrl);
}


bool keyboard::is_alt_pressed()
{
	return is_key_pressed(key_alt);
}


#endif // WIN32

GTB_END_NAMESPACE
