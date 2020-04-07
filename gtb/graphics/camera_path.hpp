
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


#ifndef GTB_CAMERA_PATH_INCLUDED
#define GTB_CAMERA_PATH_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/camera.hpp>

GTB_BEGIN_NAMESPACE


class CameraPath {
public:
	CameraPath();

	void add_camera(const Camera &c);
	void clear_cameras();
	unsigned num_cameras() const;
	const Camera &camera(unsigned i) const;
	const Camera &current_camera() const;

	void set_current_frame(unsigned i);
	void advance_frame();
	unsigned current_frame() const;

	void read(const char *file_name);
	void write(const char *file_name) const;

protected:
	int _current_frame;
	std::vector<Camera> _cameras;
};


GTB_END_NAMESPACE

#endif // GTB_CAMERA_PATH_INCLUDED
