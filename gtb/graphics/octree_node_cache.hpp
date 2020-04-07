
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


#ifndef OCTREE_NODE_CACHE_INCLUDED
#define OCTREE_NODE_CACHE_INCLUDED

#include <gtb/graphics/octree_node.hpp>
#include <list>
#include <map>
#include <string>

GTB_BEGIN_NAMESPACE


class OctreeNodeCache {
public:
	explicit OctreeNodeCache(const char *file_name_prefix);
	~OctreeNodeCache();

	// bring to memory, and assign high priority
	void fetch(OctreeNode *node);

	// bring to memory, and assign low priority
	void prefetch(OctreeNode *node, bool &node_was_in_cache);

	void reset_counters();
	unsigned num_hits() const;
	unsigned num_misses() const;
	unsigned num_prefetches() const;

	enum { DEFAULT_SIZE = 50 };
	void set_size(unsigned size);
	void double_size();
	void halve_size();
	unsigned size() const;

	const std::list<OctreeNode *> &nodes() const;
    bool contains(const std::vector<OctreeNode *> &nodes) const;

protected:
	void read_node_data(OctreeNode *node) const;

	typedef std::list<OctreeNode *> NodeList;
	typedef std::map<OctreeNode *, NodeList::iterator> NodeMap;

	std::string _file_name_prefix;
	NodeList _nodes;
	NodeMap _node_map;
	unsigned _size;
	unsigned _num_hits;
	unsigned _num_misses;
	unsigned _num_prefetches;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include "octree_node_cache.ipp"
#endif

#endif // OCTREE_NODE_CACHE_INCLUDED
