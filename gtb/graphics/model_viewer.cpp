
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
#include <gtb/graphics/model_viewer.hpp>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/model_viewer.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


ModelViewer::ModelViewer()
	: _window(-1),
	  _plp(NULL),
	  _mode(LOOK_MODE),
	  _has_picked_point(false)
{
	enable(SHOW_MODEL);
	enable(SHOW_OCTREE);
}


ModelViewer::~ModelViewer()
{
	delete _plp;
}


void ModelViewer::open_model(const char *file_name)
{
	assert(_plp == NULL);
	_plp = new Plp(file_name);
	_view.init_exterior_view(octree().bounding_box());
}


void ModelViewer::pick_vertex(int x, int y, real_t max_dist)
{
	assert(_plp != NULL);
	_world_ray = _view.world_ray(x, y);
	_has_picked_point = octree().pick_vertex(_world_ray, max_dist, _hits,
						 _picked_point, _plp->cache());
}


void ModelViewer::find_visible_set(const View &a_view)
{
	assert(_plp != NULL);
	_plp->find_visible_set(a_view);
}


void ModelViewer::render_node(const OctreeNode *node, const View &a_view) const
{
	assert(_plp != NULL);
	assert(node != NULL);
	assert(node->has_model());
	bool points_as_cones = false;
	unsigned num_points_rendered;
	node->render(a_view.camera(),
		     COLOR_RGB_WHITE,
		     COLOR_RGB_RED,
		     is_enabled(SHOW_MODEL),
		     is_enabled(SHOW_OCTREE),
		     _plp->lod(),
		     &num_points_rendered,
		     points_as_cones);
}


void ModelViewer::render_visible_set(const View &a_view) const
{
	assert(_plp != NULL);
    const std::vector<OctreeNode *> &V = _plp->visible_set();
	for (unsigned i = 0; i < V.size(); ++i) {
		OctreeNode *node = V[i];
		_plp->cache().fetch(node);
		render_node(node, a_view);
	}
}


GTB_END_NAMESPACE
