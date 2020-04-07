
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


#ifndef GTB_OCTREE_NODE_INCLUDED
#define GTB_OCTREE_NODE_INCLUDED

#include <gtb/common.hpp>
#include <gtb/graphics/coordinate_system.hpp>
#include <gtb/graphics/ray3.hpp>
#include <gtb/graphics/indexed_triangle_set.hpp>
#include <limits.h>
#include <set>

GTB_BEGIN_NAMESPACE


#define MIN_PCA_POINTS 4

class Octree;


class OctreeNode : public Box3 {
public:
	enum Direction { L, R, D, U, B, F };

	//enum Octant {LDB, RDB, LUB, RUB, LDF, RDF, LUF, RUF, ROOT };
	enum Octant {LDB, LDF, LUB, LUF, RDB, RDF, RUB, RUF, ROOT };

	OctreeNode(const Box3 &box,
		   int id,
		   Octant octant,
		   unsigned depth,
		   OctreeNode *parent);

	virtual ~OctreeNode();

	int id() const;
  	Octant octant() const;
  	unsigned depth() const;
	unsigned num_vertices() const;
	unsigned num_triangles() const;
	bool is_leaf() const;
  	OctreeNode *parent() const;

	void insert_vertices(const IndexedTriangleSet &its,
			     Octree &octree,
			     bool verbose,
			     FILE *log_file);
	void insert_triangles(const IndexedTriangleSet &its,
			      std::vector<unsigned> &replications,
			      bool verbose,
			      FILE *log_file);

	void collect_leaves(std::vector<const OctreeNode *> &leaves) const;
	void collect_leaves(std::vector<OctreeNode *> &leaves);
	void collect_leaves_excluding_subtree(
		const OctreeNode *subtree_to_exclude,
		std::vector<OctreeNode *> &leaves);
  	void collect_visible_leaves(
		const Camera &cam,
		std::vector<OctreeNode *> &nodes);
  	void collect_visible_leaves(
		const Camera &cam,
		std::vector<const OctreeNode *> &nodes) const;
	void collect_nodes(std::vector<OctreeNode *> &nodes);
	unsigned num_leaves() const;
	unsigned num_nodes() const;

	void create_submodel(const IndexedTriangleSet &its);
	void add_submodel(const IndexedTriangleSet &its);
	bool has_model() const;
	const IndexedTriangleSet &model() const;
	IndexedTriangleSet &model();
	void free_model();

	// load/save: recursive
	void load_structure(FILE *fp, Octree &octree);
	void save_structure(FILE *fp) const;

	void load_data(const char *file_name_prefix);

	// read/write: non-recursive
	void read_data(const char *file_name_prefix);
	void write_data(const char *file_name_prefix) const;
	
	void print_leaves(const char *prefix) const;
	void scan_leaves_principal_components(const char *prefix);

	Box3 bounding_box() const;
	//Point3 centroid() const;

	static unsigned max_num_vertices();
	static void set_max_num_vertices(unsigned max_num_vertices);

	static unsigned max_depth();
	static void set_max_depth(unsigned max_depth);

	static Vector3 normal(Direction d);

	OctreeNode *find_leaf_containing(const Point3 &p);
	OctreeNode *find_node(int id);

  	void insert_neighbor(OctreeNode *node, Direction dir);
	OctreeNode *neighbor(Direction direction);
	bool is_neighbor(const OctreeNode *node) const;
	void clear_neighbors();
  	const std::vector<OctreeNode *> &neighbors() const;
  	const std::vector<Direction> &neighbor_directions() const;

	int max_id() const;
	unsigned find_max_depth() const;
	bool has(OctreeNode *node) const;
	void adopt(OctreeNode *kid, Octree &octree);

	class hit {
	public:
		hit(OctreeNode *a_node, real_t a_t1, real_t a_t2);
		OctreeNode *node() const;
		real_t t1() const;
		real_t t2() const;
		friend bool operator<(const hit &h1, const hit &h2);

	protected:
		OctreeNode *m_node;
		real_t m_t1;
		real_t m_t2;
	};

	typedef std::set<hit> hit_set;

	void intersect(const Ray3 &ray, hit_set &hits);

	void reset_flags();
	void mark_as_clean();
	void mark_as_projected();
	void mark_as_enqueued();
	void mark_as_unseen();

	void init_solidity(unsigned max_num_tri);
	real_t initial_solidity() const;
	void set_solidity(real_t solidity);
	real_t current_solidity() const;
  	//friend bool operator<(const OctreeNode &a, const OctreeNode &b);
  	static bool compare_solidity(const OctreeNode *a, const OctreeNode *b);

	class CompareSolidity {
	public:
		bool operator()(const OctreeNode *a, const OctreeNode *b);
	};

	unsigned layer() const;
	void set_layer(unsigned layer);

	unsigned time_stamp() const;
	void set_time_stamp(unsigned time_stamp);

	const CoordinateSystem &coordinate_system() const;

	void init_occupancy_statistics(int &min_num_tri,
				       int &max_num_tri,
				       int &total_num_tri,
				       int &total_num_nodes) const;

	void render(const Camera &cam,
		    const ColorRgb &octree_color /*= ColorRgb::WHITE*/,
		    const ColorRgb &default_model_color /*= ColorRgb::RED*/,
		    bool show_model_enabled /*= true*/,
		    bool show_octree_enabled /*= true*/,
		    real_t lod /*= 1.0*/,
		    unsigned *num_points_rendered /*= NULL*/,
		    bool points_as_cones /*= false*/) const;
	void render_solid_box() const;
	void render_geometry() const;

	const OctreeNode *child(unsigned i) const;
	OctreeNode *child(unsigned i);

  	unsigned num_new_vertices() const;
  	unsigned num_new_triangles() const;

protected:

	enum { NUM_CHILDREN = 8 };
	static unsigned _max_num_vertices;
	static unsigned _max_depth;

	bool insert_vertex(unsigned index,
			   const std::vector<Point3> &vertices,
			   Octree &octree);
	bool insert_triangle(unsigned index,
			     const std::vector<Point3> &vertices,
			     const std::vector<IndexedTriangle> &triangles,
			     std::vector<unsigned> &replications);

	// from Box3
	void read(FILE *);
	void write(FILE *) const;

	Octant octant(const Point3 &p) const;
	Octant octant(unsigned i, unsigned j, unsigned k) const;
	Box3 child_box(unsigned i, unsigned j, unsigned k) const;
	OctreeNode *smallest_child_containing(const Box3 &box);
	void increment_depth();
	void create_children(Octree &octree);
	void insert_triangle(unsigned index, const Box3 &b,
			     std::vector<unsigned> &replications);

	// read/write: non-recursive
	void read_structure(FILE *fp);
	void write_structure(FILE *fp) const;

  	static bool adj(unsigned direction, unsigned octant);
  	static unsigned reflect(unsigned direction, unsigned octant);

	static void increment_global_ray_id();
	void rec_intersect(const Ray3 &ray, hit_set &hits);

	int _id;
	Octant _octant;
	unsigned _depth;
	unsigned _num_vertices;
	unsigned _num_triangles;
	bool _is_leaf;

	OctreeNode *_parent;
	OctreeNode *_children[NUM_CHILDREN];
	IndexedTriangleSet *_model;

	// only used when building from a model
 	std::vector<int> _vertices;
  	std::vector<int> _triangles;

  	std::vector<Direction> _neighbor_directions;
  	std::vector<OctreeNode *> _neighbors;

	real_t _initial_solidity;
	real_t _current_solidity;

	unsigned _time_stamp;
	unsigned _layer;
	CoordinateSystem _cs;

	unsigned _ray_id;
	static unsigned _global_ray_id;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/octree_node.ipp>
#endif

#endif // GTB_OCTREE_NODE_INCLUDED
