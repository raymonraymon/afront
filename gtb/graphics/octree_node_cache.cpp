
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

#include "octree_node_cache.hpp"


#ifdef OUTLINE
#define inline
#include "octree_node_cache.ipp"
#undef inline
#endif


using namespace std;

GTB_BEGIN_NAMESPACE


OctreeNodeCache::OctreeNodeCache(const char *file_name_prefix)
	: _file_name_prefix(file_name_prefix),
	  _size(DEFAULT_SIZE)
{
}


OctreeNodeCache::~OctreeNodeCache()
{
}


void OctreeNodeCache::fetch(OctreeNode *node)
{
	assert(node != NULL);
	assert(node->is_leaf());
	assert(node->num_vertices() > 0);

	// search for node in the cache
	NodeMap::iterator mi = _node_map.find(node);
	if (mi != _node_map.end()) {
		// node is in the cache
		_num_hits++;
		NodeList::iterator li = mi->second;
		assert(*li == node);
		assert(node->has_model());

		// update node's priority
		_nodes.erase(li);
		li = _nodes.insert(_nodes.end(), node);
		assert(*li == node);
		_node_map[node] = li;
	} else {
		// node is not in the cache
		_num_misses++;
		read_node_data(node);

		// replacement policy
		assert(_nodes.size() <= _size);
		if (_nodes.size() == _size) {
			//fprintf(stderr, "replacing\n");
			OctreeNode *replaced = _nodes.front();
			replaced->free_model();
			_nodes.pop_front();

			NodeMap::iterator mi2 = _node_map.find(replaced);
			assert(mi2 != _node_map.end());
			_node_map.erase(mi2);
			assert(replaced != node);
		}

		// put node in the cache
		NodeList::iterator li = _nodes.insert(_nodes.end(), node);
		assert(*li == node);
		_node_map[node] = li;
	}
}


bool OctreeNodeCache::contains(const vector<OctreeNode *> &v) const
{
	for (unsigned i = 0; i < v.size(); i++) {
		OctreeNode *node = v[i];
		NodeMap::const_iterator mi = _node_map.find(node);
		if (mi == _node_map.end()) {
			return false;
		}
	}
	return true;
}


void OctreeNodeCache::prefetch(OctreeNode *node, bool &node_was_in_cache)
{
	NodeMap::iterator mi = _node_map.find(node);
	node_was_in_cache = (mi != _node_map.end());
	if (node_was_in_cache) {
		return;
	}
	// read node's data
	read_node_data(node);

	// insert node with low priority
	NodeList::iterator li = _nodes.insert(_nodes.begin(), node);
	assert(*li == node);
	_node_map[node] = li;

	// update prefetch count (should we count hits too?)
	_num_prefetches++;
}


void OctreeNodeCache::set_size(unsigned arg_size)
{
	assert(arg_size > 0);
	if (_size != arg_size) {
		NodeList::iterator li;
		for (li = _nodes.begin(); li != _nodes.end(); li++) {
			OctreeNode *node = *li;
			node->free_model();
		}
		_nodes.clear();
		_node_map.clear();
		_size = arg_size;
	}
	assert(_nodes.size() <= _size);
}


GTB_END_NAMESPACE
