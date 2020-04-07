
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


#ifndef PLP_INCLUDED
#define PLP_INCLUDED

#include <gtb/graphics/octree.hpp>
#include <gtb/graphics/octree_node.hpp>
#include <gtb/graphics/octree_node_cache.hpp>
#include <gtb/graphics/view.hpp>
#include <gtb/graphics/image_tile_set.hpp>
#include <gtb/math/math.hpp>
#include "front.hpp"
#include "extinction_coef_set.hpp"

GTB_BEGIN_NAMESPACE


// The Prioritized-Layered Projection Algorithm for Visible Set
// Estimation

class Plp {
public:
	// visibility modes
	enum { VM_AGGRESSIVE, VM_CONSERVATIVE };

	// priority modes
	enum { PM_LAYER, PM_EXTINCTION_COEF };

	explicit Plp(const char *file_name_prefix,
		     int visibility_mode = VM_AGGRESSIVE,
		     int priority_mode = PM_LAYER);
	~Plp();

	void start(const View &view);

	// Returns true as long as the front is not empty and we
	// haven't run out of budget.
	bool step();

  	void conservative_finish();

	void find_visible_set(const View &view);

	void prefetch_extra_layers();

	const plp_front &front() const;
	const std::vector<OctreeNode *> &visible_set() const;

	const Octree &octree() const;
	Octree &octree();

	const OctreeNodeCache &cache() const;
	OctreeNodeCache &cache();

	void set_visibility_mode(int mode);
	void toggle_visibility_mode();
	int visibility_mode() const;

	enum { DEFAULT_BUDGET = 10000 };
	void set_budget(unsigned budget);
	unsigned budget() const;
	void double_budget();
	void halve_budget();

	void set_lod(float lod);
	float lod() const;
	void double_lod();
	void halve_lod();

	void set_num_extra_layers(unsigned layers);
	unsigned num_extra_layers() const;

	enum {
		NS_CLEAN,	// not visited
		NS_PROJECTED,	// rendered or culled with vis. tests
		NS_ENQUEUED,	// will be projected
		NS_UNSEEN	// out of the view frustum
	};

	typedef unsigned char node_state_t;
	node_state_t node_state(const OctreeNode *node) const;
	const char *node_state_str(const OctreeNode *node) const;

protected:
  	void front_push(OctreeNode *);
  	OctreeNode *front_pop();
	void project(OctreeNode *node);
  	bool is_potentially_visible(OctreeNode *node);
	bool reached_budget() const;
  	void add_neighbors_to_front(OctreeNode *node);
	void find_visible_front();
	float extinction_coef_priority(OctreeNode *node);
	extinction_coef get_extinction_coef(OctreeNode *node,
					    const Vector3 &v);
	void set_node_state(const OctreeNode *node, node_state_t state);
	void reset_node_states();

	plp_front _front;
	std::vector<OctreeNode *> _visible_set;

	int _visibility_mode;
	int _priority_mode;
	unsigned _budget;
	float _lod;
	unsigned _num_vertices_rendered;
	unsigned _num_triangles_rendered;
	unsigned _num_nodes_projected;
	unsigned _num_nodes_enqueued;
	unsigned _time_stamp;
	unsigned _num_extra_layers;
	ImageTileSet _tiles;

	Octree _octree;
	bool _model_is_point_set;
	OctreeNodeCache _cache;
	const View *_view;
	OctreeNode *_closest_leaf;

	class node_datum {
	public:
		node_state_t state;
		extinction_coef_set extinction_coeffs;
	};
	std::vector<node_datum> _node_data;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include "plp.ipp"
#endif

#endif // PLP_INCLUDED
