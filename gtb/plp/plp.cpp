
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

#include "plp.hpp"
#include <gtb/graphics/ogltools.h>
#include <gtb/debug/debug.hpp>
#include <gtb/error/error.hpp>


#ifdef OUTLINE
#define inline
#include "plp.ipp"
#undef inline
#endif


using namespace std;

GTB_BEGIN_NAMESPACE


Plp::Plp(const char *file_name_prefix,
	 int a_visibility_mode,
	 int a_priority_mode)
	: _visibility_mode(a_visibility_mode),
	  _priority_mode(a_priority_mode),
	  _budget(DEFAULT_BUDGET),
	  _lod(1.0),
	  _num_vertices_rendered(0),
	  _num_triangles_rendered(0),
	  _num_nodes_projected(0),
	  _num_nodes_enqueued(0),
	  _time_stamp(0),
	  _num_extra_layers(0),
	  _tiles(0, 0, 512, 512, 64, 64, GL_RGBA, GL_UNSIGNED_BYTE),
	  _cache(file_name_prefix),
	  _view(NULL),
	  _closest_leaf(NULL)
{
	_octree.load_structure(file_name_prefix);
	assert(!_octree.bounding_box().is_empty());
	_octree.init();

	if (_octree.num_triangles() == 0) {
		if (_octree.num_vertices() == 0) {
			GTB_ERROR("empty octree");
		} else {
			_model_is_point_set = true;
		}
	} else {
		_model_is_point_set = false;
	}

	_node_data.resize(_octree.num_nodes());
	if (_octree.find_max_id() >= (int) _node_data.size()) {
		GTB_ERROR("bad max id");
	}
	reset_node_states();

	if (_priority_mode == PM_EXTINCTION_COEF) {
		vector<OctreeNode *> leaves;
		_octree.collect_leaves(leaves);
		for (unsigned i = 0; i < leaves.size(); i++) {
			OctreeNode *node = leaves[i];
			assert(node != NULL);
			// FIXME: extend to point sets
			if (node->num_triangles() == 0) {
				continue;
			}
			_node_data[node->id()].extinction_coeffs.read(
				file_name_prefix, node->id());
		}
	}
}


Plp::~Plp()
{
}


void Plp::start(const View &view)
{
	_view = &view;
	assert(_view != NULL);

	_num_vertices_rendered = 0;
	_num_triangles_rendered = 0;
	_num_nodes_projected = 0;
	_num_nodes_enqueued = 0;
	_time_stamp = 0;

	reset_node_states();
	_front.clear();
	_visible_set.clear();
	_cache.reset_counters();

 	_closest_leaf = _octree.find_closest_leaf(_view->camera());
	if (_closest_leaf == NULL) {
		return;
	}
 	assert(_closest_leaf != NULL);
 	assert(_closest_leaf->is_leaf());
 	_closest_leaf->set_solidity(0.0);
 	_closest_leaf->set_layer(0);
	front_push(_closest_leaf);
}


bool Plp::step()
{
	assert(_view != NULL);

	if (_front.empty()) {
		return false;
	}
	OctreeNode *node = front_pop();
	project(node);

	// Note: we add the neighbors to the front before checking the
	// budget to guarantee that the neighbors will be on the front
	// for cPLP.
	add_neighbors_to_front(node);

	if (reached_budget()) {
		return false;
	}
	return true;
}


void Plp::front_push(OctreeNode *node)
{
	assert(node != NULL);
	assert(node->is_leaf());
	assert(node_state(node) != NS_ENQUEUED);

	unsigned n = _front.size();
	node->set_time_stamp(_time_stamp);
	_time_stamp++;
	_front.insert(node);
	if (_front.size() != (n + 1)) {
		GTB_ERROR("front_push failed: check for duplicate keys");
	}
	set_node_state(node, NS_ENQUEUED);
	_num_nodes_enqueued++;
}


// returns node most likely to be visible
OctreeNode *Plp::front_pop()
{
	OctreeNode *node = *(_front.begin());
	assert(node != NULL);
	assert(node->is_leaf());
	_front.erase(_front.begin());
	return node;
}


void Plp::project(OctreeNode *node)
{
	assert(node != NULL);
	assert(node_state(node) == NS_ENQUEUED);
	set_node_state(node, NS_PROJECTED);
	_num_nodes_projected++;
	if (node->num_vertices() == 0) {
		return;
	}
	if (is_potentially_visible(node)) {
		_visible_set.push_back(node);
		_num_vertices_rendered += (unsigned) (node->num_vertices() *
						      _lod);
		_num_triangles_rendered += node->num_triangles();
	}
}


bool Plp::is_potentially_visible(OctreeNode *node)
{
    	assert(_view != NULL);
    	assert(node != NULL);
    	return _view->camera().sees_ignoring_near_plane(node->bounding_box());
}


bool Plp::reached_budget() const
{
	if (_model_is_point_set) {
		return _num_vertices_rendered >= _budget;
	} else {
		return _num_triangles_rendered >= _budget;
	}
}


// transitive closure
void Plp::prefetch_extra_layers()
{
	assert(_view != NULL);
	assert(_cache.num_prefetches() == 0);

	if (_visible_set.size() >= _cache.size()) {
		// The cache is too small and it is thrashing, so
		// there is no point in continuing to thrash.
		return;
	}
	assert(_cache.contains(_visible_set));
	assert(_visible_set.size() < _cache.size());
	unsigned max_prefetches = _cache.size() - _visible_set.size();

	vector<OctreeNode *> v = _visible_set;
	unsigned begin = 0;
	unsigned end = v.size();
	// for each layer
	for (unsigned i = 0; i < _num_extra_layers; i++) {
		if (begin == end) {
			break;
		}
		// for each node in [begin, end)
		for (unsigned j = begin; j < end; j++) {
			assert(j < v.size());
			OctreeNode *node = v[j];

			// for each neighbor of node
			const vector<OctreeNode *> &neighbors =
				node->neighbors();
			for (unsigned k = 0; k < neighbors.size(); k++) {
				OctreeNode *nk = neighbors[k];

				// ignore empty neighbors
				if (nk->num_vertices() == 0) {
					continue;
				}

				// prefetch neighbor
				bool node_was_in_cache;
				_cache.prefetch(nk, node_was_in_cache);
				if (!node_was_in_cache) {
					v.push_back(nk);
				}

				// check for thrashing
				assert(_cache.num_prefetches()
				       <= max_prefetches);
				if (_cache.num_prefetches()
				    == max_prefetches) {
					assert(_cache.contains(_visible_set));
					return;
				}
			}
		}
		begin = end;
		end = v.size();
	}
	assert(_cache.contains(_visible_set));
}


float Plp::extinction_coef_priority(OctreeNode *node)
{
    	assert(_view != NULL);
    	assert(_closest_leaf != NULL);
	assert(node != _closest_leaf);
	assert(is_potentially_visible(node));

	Vector3 v = node->centroid() - _view->camera().origin();
	float l2 = v.squared_length();
	float k = get_extinction_coef(node, v);

	Ray3 r(_view->camera().origin(), v);
	OctreeNode::hit_set hits;
	_closest_leaf->intersect(r, hits);
	bool found = false;
	OctreeNode::hit_set::iterator hi;
	for (hi = hits.begin(); hi != hits.end(); hi++) {
		const OctreeNode::hit &hit = *hi;
		if (hit.node() == node) {
			found = true;
			break;
		}
		float k2 = get_extinction_coef(hit.node(), v);
		float travel;
		if (hit.t1() > 0.0) {
			travel = hit.t2() - hit.t1();
		} else {
			travel = hit.t2();
		}
		float dist2 = hit.t2() * hit.t2() * l2;
		k *= (1.0 - k2 * travel / dist2);
	}
	assert(found);
	return k;
}


void Plp::add_neighbors_to_front(OctreeNode *node)
{
	assert(_view != NULL);
	assert(node != NULL);

	// for each neighbor
	const vector<OctreeNode *> &neighbors = node->neighbors();
	for (unsigned i = 0; i < neighbors.size(); i++) {
		OctreeNode *neighbor = neighbors[i];

		// sanity checks
		assert(neighbor != NULL);
		assert(neighbor != node);
		assert(neighbor->is_leaf());
		//assert(neighbor->is_neighbor(node));
		//assert(node->is_neighbor(neighbor));

		// check if neighbor has been visited
		if (node_state(neighbor) != NS_CLEAN) {
			continue;
		}

		// perform view frustum culling
  		if (!is_potentially_visible(neighbor)) {
  			set_node_state(neighbor, NS_UNSEEN);
  			continue;
  		}

		// set layer
		neighbor->set_layer(node->layer() + 1);

		// set priority
		switch (_priority_mode) {
		case PM_LAYER:
			// onion-like
			neighbor->set_solidity(neighbor->layer());
			break;
		case PM_EXTINCTION_COEF:
			neighbor->set_solidity(extinction_coef_priority(
				neighbor));
			break;
		default:
			GTB_ERROR("invalid priority mode");
			break;
		}

		// enque neighbor
		front_push(neighbor);
	}
}


// item buffer approach
void Plp::find_visible_front()
{
#ifdef HP_EXT
	vector<OctreeNode *> nodes;
	while (!_front->empty()) {
  		assert(node != NULL);
		OctreeNode *node = front_pop();
		set_node_state(node, NS_CLEAN);
		nodes.push_back(node);
	}
	for (unsigned i = 0; i < nodes.size(); i++) {
		OctreeNode *node = nodes[i];
		if (node->bounding_box().is_visible()) {
			front_push(node);
		}
	}
#else  // HP_EXT

	// save masks
	GLboolean cmask[4];
	GLboolean dmask;
	glGetBooleanv(GL_COLOR_WRITEMASK, cmask);
	glGetBooleanv(GL_DEPTH_WRITEMASK, &dmask);

	// empty front, and remember nodes to render
	vector<OctreeNode *> nodes;
	while (!_front.empty()) {
		OctreeNode *node = front_pop();
  		assert(node != NULL);
		set_node_state(node, NS_CLEAN);
		nodes.push_back(node);
	}

  	// render front
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glDepthMask(GL_FALSE);
	glClear(GL_COLOR_BUFFER_BIT);
    	for (unsigned i = 0; i < nodes.size(); i++) {
  		OctreeNode *node = nodes[i];
    		assert(node != NULL);
 		// 0 is background
		unsigned j = i + 1;
		GLubyte r = j & 0x000000FF;
		GLubyte g = j & 0x0000FF00;
		GLubyte b = j & 0x00FF0000;
		glColor3ub(r, g, b);
    		node->render_solid_box();
    	}

	// determine visible front
	_tiles.read_active();
//  	fprintf(stderr, "num active tiles=%d\n",
//  		_tiles.active_tiles().size());
	list<ImageTile *>::iterator li = _tiles.active_tiles().begin();
	while (li != _tiles.active_tiles().end()) {
		ImageTile *tile = *li;
		Image &img = tile->image();
		GLubyte *p = img.pixels();
		unsigned n = img.num_pixels();
		bool tile_active = false;
		for (unsigned i = 0; i < n; i++) {
			GLubyte r = *p++;
			GLubyte b = *p++;
			GLubyte g = *p++;
			GLubyte a = *p++;
			assert(a == 255);
			(void) a;
			unsigned j = r | g << 8 | b << 16;
			if (j == 0) {
				continue;
			}
			tile_active = true;
			j--;
			if (j >= nodes.size()) {
				fprintf(stderr, "nodes size: %u\n",
					nodes.size());
				fprintf(stderr, "bogus j: %u\n", j);
			}
			assert(j < nodes.size());
			OctreeNode *node = nodes[j];
			if (node_state(node) != NS_ENQUEUED) {
				front_push(node);
			}
		}
		list<ImageTile *>::iterator t = li;
		li++;
		if (!tile_active) {
			_tiles.deactivate(t);
		}
	}

	// restore masks
	glColorMask(cmask[0], cmask[1], cmask[2], cmask[3]);
	glDepthMask(dmask);
#endif // HP_EXT
}


void Plp::conservative_finish()
{
	assert(_view != NULL);
	const Viewport &v = _view->viewport();
	_tiles.resize(v.x_min(), v.y_min(), v.width(), v.height());
	_tiles.activate_all();

	// save state
  	glPushAttrib(GL_ENABLE_BIT);
	GLboolean cmask[4];
	GLboolean dmask;
	glGetBooleanv(GL_COLOR_WRITEMASK, cmask);
	glGetBooleanv(GL_DEPTH_WRITEMASK, &dmask);

  	glDisable(GL_LIGHTING);

        // compute the z-buffer for the approximate visible set
        glClear(GL_DEPTH_BUFFER_BIT);
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        glDepthMask(GL_TRUE);
        for (unsigned i = 0; i < _visible_set.size(); i++) {
                OctreeNode *node = _visible_set[i];
		assert(node->num_vertices() > 0);
                _cache.fetch(node);
                assert(node->has_model());
                node->render_geometry();
        }

	while (!_front.empty()) {

		// find visible front
		find_visible_front();

		// empty front, and remember nodes to project
		vector<OctreeNode *> nodes_to_project;
		while (!_front.empty()) {
			OctreeNode *node = front_pop();
			nodes_to_project.push_back(node);
		}

		// project nodes, and determine new front
		for (unsigned i = 0; i < nodes_to_project.size(); i++) {
			OctreeNode *node = nodes_to_project[i];

			project(node);
			add_neighbors_to_front(node);

			if (node->num_vertices() == 0) {
				continue;
			}

			_cache.fetch(node);
			assert(node->has_model());
			node->render_geometry();
		}
	}

	// restore state
  	glPopAttrib();
	glColorMask(cmask[0], cmask[1], cmask[2], cmask[3]);
	glDepthMask(dmask);

	print_gl_errors();
}


void Plp::find_visible_set(const View &view)
{
	start(view);
	while (step()) {
		;		// loop
	}
	if (_visibility_mode == VM_CONSERVATIVE) {
		conservative_finish();
	}
}


GTB_END_NAMESPACE
