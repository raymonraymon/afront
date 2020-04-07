
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


#ifndef GTB_MODEL_VIEWER_INCLUDED
#define GTB_MODEL_VIEWER_INCLUDED

#include <gtb/common.hpp>
#include <gtb/ui/ui.hpp>
#include <gtb/graphics/view.hpp>
#include <gtb/graphics/octree.hpp>
#include <gtb/plp/plp.hpp>
#include <map>


GTB_BEGIN_NAMESPACE


class ModelViewer {
public:
	enum Capability { SHOW_MODEL, SHOW_OCTREE, NUM_CAPS };

	enum Mode { LOOK_MODE, PICK_MODE, NUM_MODES };

	ModelViewer();
	virtual ~ModelViewer();

	void open_model(const char *file_name);

	void set_window(int win);
	int window() const;

	Mouse &mouse();
	const Mouse &mouse() const;

	View &view();
	const View &view() const;

	Octree &octree();
	const Octree &octree() const;

	bool has_model() const;

	Plp &plp();
	const Plp &plp() const;

	void find_visible_set(const View &view);
	void render_visible_set(const View &view) const;

	void enable(Capability cap);
	void disable(Capability cap);
	void toggle(Capability cap);
	bool is_enabled(Capability cap) const;

	void add_point(const Point3 &p);
	unsigned num_points() const;
	const Point3 &point(unsigned i) const;
	Point3 &point(unsigned i);
	const std::vector<Point3> &points() const;
	std::vector<Point3> &points();

	void set_mode(unsigned mode);
	unsigned mode() const;

	void pick_vertex(int x, int y, real_t max_dist);
	bool has_picked_point() const;
	const Point3 &picked_point() const;
	const Ray3 &world_ray() const;
	const OctreeNode::hit_set &hits() const;

	void post_redisplay() const;

protected:
	void render_node(const OctreeNode *node, const View &view) const;

	int _window;
	Mouse _mouse;
	View _view;
	Plp *_plp;
	bool _caps[NUM_CAPS];
	std::vector<Point3> _points;
	int _mode;
	bool _has_picked_point;
	Point3 _picked_point;
	Ray3 _world_ray;
	OctreeNode::hit_set _hits;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/model_viewer.ipp>
#endif

#endif // GTB_MODEL_VIEWER_INCLUDED
