
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
#include <gtb/graphics/image_viewer.hpp>
#include <gtb/graphics/text.hpp>
#include <gtb/error/error.hpp>
#endif // WIN32

using namespace std;

#ifdef OUTLINE
#define inline
#include <gtb/graphics/image_viewer.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


ImageViewer::ImageViewer()
	: _image(NULL),
	  _window(-1),
	  _mode(LOOK_MODE),
	  _zoom_factor(1.0),
	  _picked_point_index(-1),
	  _pick_radius(3)
{
	enable(SHOW_IMAGE);
}


ImageViewer::~ImageViewer()
{
}


void ImageViewer::set_image(const Image &img)
{
	_image = &img;
	_origin.reset(0, 0);
	_view.init_image_view(_image->width(), _image->height());
}


void ImageViewer::pick(int x, int y)
{
	_picked_point_index = -1;
	int vp_w = _view.viewport().width();
	int vp_h = _view.viewport().height();
	if ((x >= 0) && (x < vp_w)
	    && (y >= 0) && (y < vp_h)) {
		Point2 vp_p(x, y);
		Point2 img_p = viewport_to_image(vp_p);
		for (unsigned i = 0; i < num_points(); i++) {
			const Point2 &p = point(i);
			if ((fabs(p.x() - img_p.x()) <= _pick_radius)
			    && (fabs(p.y() - img_p.y()) <= _pick_radius)) {
				_picked_point_index = i;
				break;
			}
		}
	}
}


void ImageViewer::draw_image() const
{
	assert(_image != NULL);
	const Viewport &vp = _view.viewport();
	int w = min(vp.width(), (int) _image->width() - (int) _origin.x());
	int h = min(vp.height(), (int) _image->height() - (int) _origin.y());
	glPushAttrib(GL_PIXEL_MODE_BIT);
	glPixelZoom(_zoom_factor, _zoom_factor);
	_image->draw((int) _origin.x(), (int) _origin.y(), w, h);
	glPopAttrib();
}


void ImageViewer::draw_points() const
{
	float r = _pick_radius * _zoom_factor;
	for (unsigned i = 0; i < num_points(); i++) {
		const Point2 &p = point(i);
		float x = (p.x() - _origin.x()) * _zoom_factor;
		float y = (p.y() - _origin.y()) * _zoom_factor;
		if ((int) i == _picked_point_index) {
			COLOR_RGB_RED.load();
		} else {
			COLOR_RGB_GREEN.load();
		}

		// draw cross at feature
		glBegin(GL_LINES);
		glVertex2f(x - 2 * r, y);
		glVertex2f(x + 2 * r, y);
		glVertex2f(x, y - 2 * r);
		glVertex2f(x, y + 2 * r);
		glEnd();

		// draw square around feature
		glBegin(GL_LINE_LOOP);
		glVertex2f(x - r, y - r);
		glVertex2f(x + r, y - r);
		glVertex2f(x + r, y + r);
		glVertex2f(x - r, y + r);
		glEnd();
	}
}


void ImageViewer::draw_point_indices() const
{
	const Viewport &vp = _view.viewport();
	draw_text_begin(vp.x_min(), vp.width(), vp.y_min(), vp.height());
	for (unsigned i = 0; i < num_points(); i++) {
		const Point2 &p = point(i);
		if ((int) i == _picked_point_index) {
			COLOR_RGB_RED.load();
		} else {
			COLOR_RGB_GREEN.load();
		}
		int x = (int) ((p.x() - _origin.x() + 10) * _zoom_factor);
		int y = (int) ((p.y() - _origin.y() + 10) * _zoom_factor);
		draw_int(i, x, y);
	}
	draw_text_end();
}


void ImageViewer::on_mouse(int button, int state, int x, int y)
{
	bool redraw = false;
	_mouse.set_button_state(button, state);
	_mouse.set_position(x, y);

	switch (_mode) {
	case LOOK_MODE:
		break;
	case PICK_MODE:
		if ((button == GLUT_LEFT_BUTTON)
		    && (state == GLUT_DOWN)) {
			pick(x, y);
			redraw = true;
		}
		break;
	default:
		GTB_ERROR("invalid mode");
		break;
	}
	if (redraw) {
		post_redisplay();
	}
}


void ImageViewer::on_motion(int x, int y)
{
	bool redraw = false;
	_mouse.set_position(x, y);
	switch (_mode) {
	case LOOK_MODE:
		assert(_image != NULL);
		_origin.translate(-_mouse.dx() / _zoom_factor,
				  _mouse.dy() / _zoom_factor);
		clamp_origin();
		redraw = true;
		break;
	case PICK_MODE:
		assert(_image != NULL);
		if (has_picked_point()) {
			picked_point().translate(_mouse.dx() / _zoom_factor,
						 _mouse.dy() / _zoom_factor);
			clamp_point(picked_point());
		}
		redraw = true;
		break;
	default:
		GTB_ERROR("invalid mode");
		break;
	}
	if (redraw) {
		post_redisplay();
	}
}


void ImageViewer::load_points(const char *file_name)
{
	FILE *fp = fopen(file_name, "rb");
	REQUIRE(fp != NULL);
	float x, y;
	while (fscanf(fp, "%f %f", &x, &y) == 2) {
		add_point(Point2(x, y));
	}
	fclose(fp);
}


void ImageViewer::save_points(const char *file_name) const
{
	FILE *fp = fopen(file_name, "wb");
	REQUIRE(fp != NULL);
	for (unsigned i = 0; i < num_points(); i++) {
		const Point2 &p = point(i);
		fprintf(fp, "%g %g\n", p.x(), p.y());
	}
	fclose(fp);
}


void ImageViewer::go_to_next_point()
{
	if (_picked_point_index == -1) {
		_picked_point_index = 0;
	} else {
		_picked_point_index = (_picked_point_index + 1) % num_points();
	}
	center_on(picked_point());
}


void ImageViewer::go_to_previous_point()
{
	if (_picked_point_index == -1) {
		_picked_point_index = 0;
	} else {
		_picked_point_index = (_picked_point_index + num_points() - 1)
			% num_points();
	}
	center_on(picked_point());
}


void ImageViewer::center_on(const Point2 &p)
{
	assert(_image != NULL);
	const Viewport &vp = _view.viewport();
	_origin.reset(p.x() - (vp.width() / 2) / _zoom_factor,
		      p.y() - (vp.height() / 2) / _zoom_factor);
	clamp_origin();
}


void ImageViewer::clamp_origin()
{
	assert(_image != NULL); 
	const Viewport &vp = _view.viewport();
	float x = clamp((float) _origin.x(),
			(float) 0.0,
			(float) (_image->width() - 1 - 
				 vp.width() / _zoom_factor));
	float y = clamp((float) _origin.y(),
			(float) 0.0,
			(float) (_image->height() - 1 -
				 vp.height() / _zoom_factor));
	_origin.reset(x, y);
}


void ImageViewer::clamp_point(Point2 &p) const
{
	assert(_image != NULL);
	float x = clamp((float) p.x(),
			(float) 0.0,
			(float) _image->width());
	float y = clamp((float) p.y(),
			(float) 0.0,
			(float) _image->height());
	p.reset(x, y);
}


GTB_END_NAMESPACE
