
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


#ifndef GTB_OCTREE_INCLUDED
#define GTB_OCTREE_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/octree_node.hpp>
#include <gtb/graphics/octree_node_cache.hpp>
#include <gtb/graphics/indexed_triangle_set.hpp>
#include <string>
#include <map>


GTB_BEGIN_NAMESPACE


class Octree {
public:
	Octree();
	virtual ~Octree();

	void build(const IndexedTriangleSet &model,
		   FILE *log_file);
	void build(const IndexedTriangleSet &model,
		   const Box3 &bbox,
		   FILE *log_file);
	void add_model(const IndexedTriangleSet &model);

	// add_model_out_of_core() works for point sets and indexed
	// triangle sets.
	//
	// You: So why does it take an indexed triangle set?
	//
	// Me: It's because there is no class for a point set.  A point
	// set is an indexed triangle set without triangles.
	//
	// You: This is ugly!
	//
	// Me: I know...
	void add_model_out_of_core(const IndexedTriangleSet &model,
				   const char *file_name_prefix,
				   FILE *log_file);

	void load_structure(const char *file_name_prefix);
	void save_structure(const char *file_name_prefix,
			    FILE *log_file) const;

	void load_data(const char *file_name_prefix);
	void save_data(const char *file_name_prefix,
		       const IndexedTriangleSet &model,
		       FILE *log_file) const;

	void print_leaves(const char *prefix) const;
	void scan_leaves_principal_components(const char *prefix);

	void collect_leaves(std::vector<const OctreeNode *> &leaves) const;
	void collect_leaves(std::vector<OctreeNode *> &leaves) const;
	void collect_visible_leaves(
		const Camera &cam,
		std::vector<OctreeNode *> &nodes);
	void collect_visible_leaves(
		const Camera &cam,
		std::vector<const OctreeNode *> &nodes) const;
	void collect_nodes(std::vector<OctreeNode *> &nodes);
	unsigned num_leaves() const;
	unsigned num_nodes() const;

	const OctreeNode &root() const;
	Box3 bounding_box() const;
	Point3 centroid() const;
	bool contains(const Point3 &p) const;
	bool contains_all_vertices(const IndexedTriangleSet &model) const;
	OctreeNode *find_leaf_containing(const Point3 &p);
	OctreeNode *find_closest_leaf(const Camera &cam);
	OctreeNode *find_node(int id);

	void init();

	void render(const Camera &cam,
		    const ColorRgb &octree_color /*= ColorRgb::WHITE*/,
		    const ColorRgb &default_model_color /*= ColorRgb::RED*/,
		    bool show_model_enabled /*= true*/,
		    bool show_octree_enabled /*= true*/,
		    real_t lod /*= 1.0*/,
		    bool points_as_cones /*= false*/) const;

	static unsigned num_triangles_rendered();

	bool has(OctreeNode *node) const;

	void intersect(const Ray3 &ray, OctreeNode::hit_set &hits);

	bool pick(const Ray3 &ray,
		  const char *model_name,
		  OctreeNode::hit_set &hits,
		  Point3 &picked_point,
		  OctreeNode *&picked_point_node);

	bool pick_vertex(const Ray3 &ray,
			 real_t max_dist,
			 OctreeNode::hit_set &hits,
			 Point3 &picked_vertex);

	bool pick_vertex(const Ray3 &ray,
			 real_t max_dist,
			 OctreeNode::hit_set &hits,
			 Point3 &picked_vertex,
			 OctreeNodeCache &cache);

	int make_id();

	static std::string get_directory_name(const char *file_name_prefix);

	static std::string get_structure_file_name(
		const char *file_name_prefix);

	static std::string get_data_file_name(const char *file_name_prefix,
					      int id);

	static FILE *open_structure_file(const char *file_name_prefix,
					 const char *mode);

	static FILE *open_data_file(const char *file_name_prefix,
				    int id,
				    const char *mode);

	static IndexedTriangleSet *read_data(const char *file_name_prefix,
					     int id);

	static void write_data(const IndexedTriangleSet &model,
			       const char *file_name_prefix,
			       int id);

	unsigned num_vertices() const;
	unsigned num_triangles() const;
	unsigned find_max_depth() const;
	int find_max_id() const;

	unsigned find_max_num_triangles_per_leaf() const;
	float find_avg_num_triangles_per_leaf() const;

protected:
	// PLP functions
	void init_occupancy_statistics();
	void init_solidity();

	// out-of-core functions
	void grow_to_contain(const IndexedTriangleSet &model);
	void grow_to_contain(const Point3 &p);
	void commit_out_of_core_insertions(
		const std::vector<OctreeNode *> &old_leaves,
		const char *file_name_prefix,
		const IndexedTriangleSet &new_model);
	void commit_out_of_core_insertions(
		OctreeNode *old_leaf,
		const char *file_name_prefix,
		const IndexedTriangleSet &new_model);

	// other functions
	void compute_arrangement();
  	void classify_point(const Point3 &p,
			    OctreeNode::Direction dirs[3],
			    int &num_dirs);

//  	const std::vector<OctreeNode *> &outside_nodes(
//  		OctreeNode::Direction i) const;

	OctreeNode *_root;
	std::vector<OctreeNode *> _outside_nodes[6];
	unsigned _max_num_vertices;
	unsigned _max_num_triangles;
	int _max_id;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/octree.ipp>
#endif

#endif // GTB_OCTREE_INCLUDED
