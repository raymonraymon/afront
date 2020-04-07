
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
#include <gtb/graphics/octree.hpp>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/polygon3.hpp>
#include <gtb/error/error.hpp>
#include <gtb/io/io.hpp>
#include <gtb/util/util.hpp>
#include <gtb/graphics/ogltools.h>
#include <sstream>
#include <unistd.h>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/octree.ipp>
#undef inline
#endif


using namespace std;


GTB_BEGIN_NAMESPACE


// --- constructor and destructor ---

Octree::Octree()
	: _root(NULL),
	  _max_num_vertices(0),
	  _max_num_triangles(0),
	  _max_id(-1)
{
}


Octree::~Octree()
{
	if (_root != NULL) {
		delete _root;
	}
}


// --- IO functions ---

string Octree::get_directory_name(const char *file_name_prefix)
{
	ostringstream ost;
	ost << file_name_prefix << ".d";
	return ost.str();
}


string Octree::get_structure_file_name(const char *file_name_prefix)
{
	ostringstream ost;
	ost << file_name_prefix << ".d/oct";
	return ost.str();
}


string Octree::get_data_file_name(const char *file_name_prefix, int id)
{
	ostringstream ost;
	ost << file_name_prefix << ".d/" << id << ".dat";
	return ost.str();
}


FILE *Octree::open_structure_file(const char *file_name_prefix,
				  const char *mode)
{
	string file_name = get_structure_file_name(file_name_prefix);
	FILE *fp = xfopen(file_name.c_str(), mode);
	//fprintf(stderr, "opened structure file\n");
	return fp;
}


FILE *Octree::open_data_file(const char *file_name_prefix,
			     int id,
			     const char *mode)
{
	string file_name = get_data_file_name(file_name_prefix, id);
	FILE *fp = xfopen(file_name.c_str(), mode);
	//fprintf(stderr, "opened data file %s\n", file_name.c_str());
	return fp;
}


IndexedTriangleSet *Octree::read_data(const char *file_name_prefix, int id)
{
	FILE *fp = open_data_file(file_name_prefix, id, "rb");
	IndexedTriangleSet *model = new IndexedTriangleSet();
	model->read_binary(fp);
	fclose(fp);
	return model;
}


void Octree::write_data(const IndexedTriangleSet &model,
			const char *file_name_prefix,
			int id)
{
	FILE *fp = open_data_file(file_name_prefix, id, "wb");
	model.write_binary(fp);
	fclose(fp);
}


void Octree::load_structure(const char *file_name_prefix)
{
	fprintf(stderr, "reading octree structure...");
	assert(_root == NULL);
	_root = new OctreeNode(
		Box3(),		// box
		make_id(),	// id
		OctreeNode::ROOT, // octant
		0,		// depth
		NULL);		// parent

	FILE *fp = open_structure_file(file_name_prefix, "rb");
	read_unsigned(&_max_num_vertices, fp);
	read_unsigned(&_max_num_triangles, fp);
	read_int(&_max_id, fp);

	assert(_root != NULL);
	_root->load_structure(fp, *this);
	compute_arrangement();

	assert(_max_id == find_max_id());
	fclose(fp);
	fprintf(stderr, "\n");
}


void Octree::save_structure(const char *file_name_prefix,
			    FILE *log_file) const
{
	timer t;
	t.start();

	fprintf(stderr, "writing octree structure...");
	FILE *fp = open_structure_file(file_name_prefix, "wb");
	write_unsigned(_max_num_vertices, fp);
	write_unsigned(_max_num_triangles, fp);
	write_int(_max_id, fp);

	assert(_root != NULL);
	_root->save_structure(fp);

	fclose(fp);
	fprintf(stderr, "\n");

	t.stop();
	if (log_file != NULL) {
		fprintf(log_file, "save structure: %.3fs\n",
			t.elapsed_seconds());
	}
}


void Octree::load_data(const char *file_name_prefix)
{
	fprintf(stderr, "reading octree data...");
	assert(_root != NULL);
	_root->load_data(file_name_prefix);
	fprintf(stderr, "\n");
}


void Octree::save_data(const char *file_name_prefix,
		       const IndexedTriangleSet &model,
		       FILE *log_file) const
{
	timer t;
	t.start();

	fprintf(stderr, "writing octree data...");
	vector<OctreeNode *> leaves;
	collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		assert(leaf->is_leaf());
		assert(!leaf->has_model());
		if ((leaf->num_new_triangles() > 0) ||
		    (leaf->num_new_vertices() > 0)) {
			leaf->create_submodel(model);
			leaf->write_data(file_name_prefix);
			leaf->free_model();
		}
		if (((i + 1) % 100) == 0) {
			fprintf(stderr, "\rwriting octree data... %d/%d",
				i + 1, leaves.size());
		}
	}
	fprintf(stderr, "\rwriting octree data... %d/%d\n",
		leaves.size(), leaves.size());

	t.stop();
	if (log_file != NULL) {
		fprintf(log_file, "save data: %.3fs\n",
			t.elapsed_seconds());
	}
}


// --- PLP functions ---

void Octree::init_occupancy_statistics()
{
	int min_num_tri = INT_MAX;
	int max_num_tri = -1;
	int total_num_tri = 0;
	int total_num_nodes = 0;
	_root->init_occupancy_statistics(min_num_tri, max_num_tri,
					 total_num_tri, total_num_nodes);
	_max_num_triangles = max_num_tri;
}


void Octree::init_solidity()
{
	assert(_root != NULL);
	_root->init_solidity(_max_num_triangles);
}


// --- in-core building/adding ---

void Octree::build(const IndexedTriangleSet &model,
		   FILE *log_file)
{
	build(model, model.bounding_box(), log_file);
}


void Octree::build(const IndexedTriangleSet &model,
		   const Box3 &bbox,
		   FILE *log_file)
{
	_root = new OctreeNode(
		bbox,		// box
		make_id(),	// id
		OctreeNode::ROOT, // octant
		0,		// depth
		NULL);		// parent
	assert(_root != NULL);

	_root->insert_vertices(model, *this, true, log_file);

	vector<unsigned> replications;
	_root->insert_triangles(model, replications, true, log_file);

	bool show_replication_enabled = false;
	if (show_replication_enabled && (model.num_triangles() > 0)) {
		// find maximum number of replications
		unsigned max_replications = 0;
		for (unsigned i = 0; i < replications.size(); i++) {
			if (replications[i] > max_replications) {
				max_replications = replications[i];
			}
		}

		// compute histogram of replications
		vector<unsigned> histogram(max_replications + 1, 0);
		for (unsigned i = 0; i < replications.size(); i++) {
			histogram[replications[i]]++;
		}
	
		// print histogram of replications
		fprintf(stderr, "replications:\n");
		unsigned n = model.num_triangles();
		for (unsigned i = 0; i < histogram.size(); i++) {
			fprintf(stderr, "%3d %10d %9.5f%%\n", i, histogram[i],
				(float) histogram[i] / (float) n * 100.0);
		}

		// print average replication
		unsigned sum = 0;
		for (unsigned i = 0; i < histogram.size(); i++) {
			sum += i * histogram[i];
		}
		fprintf(stderr, "average replication: %f\n",
			(float) sum / (float) n);
	}
	
	//fprintf(stderr, "_max_id after insertions=%d\n", _max_id);

	compute_arrangement();
}


void Octree::add_model(const IndexedTriangleSet &model)
{
	assert(_root != NULL);

	fprintf(stderr, "adding model:\n");
	fprintf(stderr, "growing to contain model...");
	grow_to_contain(model);
	fprintf(stderr, "\n");

	_root->insert_vertices(model, *this, false, NULL);
	vector<unsigned> replications;
	_root->insert_triangles(model, replications, false, NULL);

	fprintf(stderr, "updating submodels...");
	vector<OctreeNode *> leaves;
	_root->collect_leaves(leaves);
	assert(leaves.size() > 0);
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		if ((leaf->num_new_triangles() == 0)
		    && (leaf->num_new_vertices() == 0)) {
			continue;
		}
		if (!leaf->has_model()) {
			leaf->create_submodel(model);
		} else {
			leaf->add_submodel(model);
		}
		if (((i + 1) % 1000) == 0) {
			fprintf(stderr, "\rupdating submodels... %d/%d",
				i + 1, leaves.size());
		}
	}
	fprintf(stderr, "\rupdating submodels... %d/%d\n",
		leaves.size(), leaves.size());

	compute_arrangement();
}


// --- out-of-core adding ---

void Octree::add_model_out_of_core(const IndexedTriangleSet &model,
				   const char *file_name_prefix,
				   FILE *log_file)
{
	fprintf(stderr, "adding model out of core:");
	assert(_root != NULL);

	// save old root
	OctreeNode *old_root = _root;

	// colect old leaves
  	vector<OctreeNode *> old_leaves;
  	collect_leaves(old_leaves);

	// sanity check
	for (unsigned i = 0; i < old_leaves.size(); i++) {
		OctreeNode *leaf = old_leaves[i];
		assert(!leaf->has_model());
		(void) leaf;
	}

	// make sure model fits
	grow_to_contain(model);

	// insert vertices and triangles
	_root->insert_vertices(model, *this, true, log_file);
	vector<unsigned> replications;
	_root->insert_triangles(model, replications, true, log_file);

	// sanity check
	vector<OctreeNode *> leaves;
	collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(!leaf->has_model());
		(void) leaf;
	}

// save data

	timer t;
	t.start();

	// commit insertions in the old leaves
	commit_out_of_core_insertions(old_leaves, file_name_prefix, model);

	// save new nodes created while growing
	leaves.clear();
	_root->collect_leaves_excluding_subtree(old_root, leaves);
	fprintf(stderr, "saving new leaves...");
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		if ((leaf->num_new_triangles() == 0)
		    && (leaf->num_new_vertices() == 0)) {
			continue;
		}
		assert(!leaf->has_model());
		leaf->create_submodel(model);
		leaf->write_data(file_name_prefix);
		leaf->free_model();
		if (((i + 1) % 100) == 0) {
			fprintf(stderr, "\rsaving new leaves... %d/%d",
				i + 1, leaves.size());
		}
	}
	fprintf(stderr, "\rsaving new leaves... %d/%d\n",
		leaves.size(), leaves.size());

	compute_arrangement();

	// check for memory leaks
	vector<OctreeNode *> nodes;
	collect_nodes(nodes);
	for (unsigned i = 0; i < nodes.size(); i++) {
		assert(!nodes[i]->has_model());
	}

	t.stop();
	if (log_file != NULL) {
		fprintf(log_file, "save data: %.3fs\n",
			t.elapsed_seconds());
	}
}


//  static void print_box(const debug &d, const char *msg, const Box3 &b)
//  {
//  	d.print("%s: (%f, %f, %f) (%f, %f, %f)\n", msg,
//  		b.x_min(), b.y_min(), b.z_min(),
//  		b.x_max(), b.y_max(), b.z_max());
//  }



bool Octree::contains_all_vertices(const IndexedTriangleSet &model) const
{
	for (unsigned i = 0; i < model.num_vertices(); i++) {
		if (!contains(model.vertex(i))) {
			return false;
		}
	}
	return true;
}


void Octree::grow_to_contain(const IndexedTriangleSet &model)
{
	// get bbox
	const Box3 &b = model.bounding_box();

	// min point
	grow_to_contain(b.min_point());

	// max point
	grow_to_contain(b.max_point());

	// sanity check
	assert(contains_all_vertices(model));
}


void Octree::grow_to_contain(const Point3 &p)
{
    const Vector3 &x = Vector3::VECTOR3_POSITIVE_X;
	const Vector3 &y = Vector3::VECTOR3_POSITIVE_Y;
	const Vector3 &z = Vector3::VECTOR3_POSITIVE_Z;

	while (!contains(p)) {
		// get current dimensions
		real_t lx = _root->x_length();
		real_t ly = _root->y_length();
		real_t lz = _root->z_length();
		
		// get vector from p to current centroid
		Point3 c = centroid();
		Vector3 v = p - c;

		// get signs of projections of v on the axes
		real_t sx = (v.dot(x) >= 0.0) ? 1.0 : -1.0;
		real_t sy = (v.dot(y) >= 0.0) ? 1.0 : -1.0;
		real_t sz = (v.dot(z) >= 0.0) ? 1.0 : -1.0;

		// get centroid of new root
		Vector3 dx = ((real_t)0.5) * sx * lx * x;
		Vector3 dy = ((real_t)0.5) * sy * ly * y;
		Vector3 dz = ((real_t)0.5) * sz * lz * z;
		Point3 nrc = c + dx + dy + dz;

		// create new root
		Box3 nrb(nrc, 2.0 * lx, 2.0 * ly, 2.0 * lz);
		OctreeNode *nr = new OctreeNode(
			nrb,	// box
			make_id(), // id
			OctreeNode::ROOT, // octant
			0,	// depth;
			NULL);	// parent
		nr->adopt(_root, *this);
		_root = nr;

		assert(!_root->is_leaf());
		for (unsigned i = 0; i < 8; i++) {
			assert(_root->child(i) != NULL);
		}
	}
}


void Octree::commit_out_of_core_insertions(
	const vector<OctreeNode *> &old_leaves,
	const char *file_name_prefix,
	const IndexedTriangleSet &new_model)
{
	// for each old leaf
	fprintf(stderr, "saving old leaves...");
  	for (unsigned i = 0; i < old_leaves.size(); i++) {
  		OctreeNode *old_leaf = old_leaves[i];

		// sanity check
		assert(!old_leaf->has_model());

		// trivial case: no insertion was made
		if (old_leaf->is_leaf()
		    && (old_leaf->num_new_vertices() == 0)
		    && (old_leaf->num_new_triangles() == 0)) {
			continue;
		}

		commit_out_of_core_insertions(old_leaf,
					      file_name_prefix,
					      new_model);

		if (((i + 1) % 100) == 0) {
			fprintf(stderr, "\rsaving old leaves... %d/%d",
				i + 1, old_leaves.size());
		}

 	}
	fprintf(stderr, "\rsaving old leaves... %d/%d\n",
		old_leaves.size(), old_leaves.size());
}


// old_leaf may not be a leaf anymore
// old data is on disk
// new leaves contain indices to new data being added
void Octree::commit_out_of_core_insertions(
	OctreeNode *old_leaf,
	const char *file_name_prefix,
	const IndexedTriangleSet &new_model)
{
	// get old data
	IndexedTriangleSet *tmp_model;
	if (old_leaf->num_vertices() > 0) {
		tmp_model = Octree::read_data(file_name_prefix,
					      old_leaf->id());
	} else {
		tmp_model = new IndexedTriangleSet;
	}

	// collect leaves
	vector<OctreeNode *> leaves;
	old_leaf->collect_leaves(leaves);
	assert(leaves.size() > 0);

	// merge old and new data
	unsigned nt = tmp_model->num_triangles();
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		if ((leaf->num_new_triangles() == 0)
		    && (leaf->num_new_vertices() == 0)) {
			continue;
		}
		nt += leaf->num_new_triangles();
		leaf->create_submodel(new_model);
		tmp_model->add(leaf->model());
		leaf->free_model();
	}
	assert(tmp_model->is_valid());
	assert(tmp_model->num_triangles() == nt);

	// reinsert everything
	old_leaf->insert_vertices(*tmp_model, *this, false, NULL);
	vector<unsigned> replications;
	old_leaf->insert_triangles(*tmp_model, replications, false, NULL);

	// recollect leaves
	leaves.clear();
	old_leaf->collect_leaves(leaves);

	// create submodels
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		if ((leaf->num_new_triangles() == 0)
		    && (leaf->num_new_vertices() == 0)) {
			continue;
		}
		assert(!leaf->has_model());
		leaf->create_submodel(*tmp_model);
		leaf->write_data(file_name_prefix);
		leaf->free_model();
	}

	// remove old data file, if necessary
	if (!old_leaf->is_leaf()) {
		string file_name = get_data_file_name(file_name_prefix,
						      old_leaf->id());
		unlink(file_name.c_str());
	}

	// clean up
	delete tmp_model;
	tmp_model = NULL;
}


// --- other funtions ---

inline bool depth_sort_node(const OctreeNode *a, const OctreeNode *b)
{
	assert(a != NULL);
	assert(b != NULL);
	return a->depth() < b->depth();
}


void Octree::compute_arrangement()
{
	// clear outside nodes
	for (unsigned i = 0; i < 6; i++) {
		_outside_nodes[i].clear();
	}

	// clear neighbors
	vector<OctreeNode *> nodes;
	collect_nodes(nodes);
	for (unsigned i = 0; i < nodes.size(); i++) {
		assert(nodes[i] != NULL);
		nodes[i]->clear_neighbors();
	}

	// sort leaves by depth
	nodes.clear();
	collect_leaves(nodes);
	sort(nodes.begin(), nodes.end(), depth_sort_node);
	for (unsigned i = 0; i < nodes.size() - 1; i++) {
		assert(nodes[i]->depth() <= nodes[i + 1]->depth());
	}

	// for each leaf
	for (unsigned i = 0; i < nodes.size(); i++) {
		OctreeNode *node = nodes[i];
		assert(node != NULL);
		assert(node->is_leaf());

		static OctreeNode::Direction dir[6] = {
			OctreeNode::L,
			OctreeNode::R,
			OctreeNode::D,
			OctreeNode::U,
			OctreeNode::B,
			OctreeNode::F
		};
		static OctreeNode::Direction opp_dir[6] = {
			OctreeNode::R,
			OctreeNode::L,
			OctreeNode::U,
			OctreeNode::D,
			OctreeNode::F,
			OctreeNode::B
		};

		// for each direction
		for (unsigned d = 0; d < 6; d++) {
			// get neighbor
			OctreeNode *x = node->neighbor(dir[d]);
			if (x == NULL) {
				// no neighbor: external node
				_outside_nodes[dir[d]].push_back(node);
			} else if (x->is_leaf()) {
				if (!node->is_neighbor(x)) {
					node->insert_neighbor(x, dir[d]);
				}
				if (!x->is_neighbor(node)) {
					x->insert_neighbor(node, opp_dir[d]);
				}
			}
		}
	}
}


void Octree::classify_point(const Point3 &p,
			    OctreeNode::Direction dirs[3],
			    int &num_dirs)
{
	num_dirs = 0;
	const Point3 &bb_min = _root->min_point();
	const Point3 &bb_max = _root->max_point();

	if (p[0] < bb_min[0]) {
		dirs[num_dirs] = OctreeNode::L;
		++num_dirs;
	} else if (p[0] > bb_max[0]) {
		dirs[num_dirs] = OctreeNode::R;
		++num_dirs;
	}

	if (p[1] < bb_min[1]) {
		dirs[num_dirs] = OctreeNode::D;
		++num_dirs;
	} else if (p[1] > bb_max[1]) {
		dirs[num_dirs] = OctreeNode::U;
		++num_dirs;
	}

	if (p[2] < bb_min[2]) {
		dirs[num_dirs] = OctreeNode::B;
		++num_dirs;
	} else if (p[2] > bb_max[2]) {
		dirs[num_dirs] = OctreeNode::F;
		++num_dirs;
	}
}


OctreeNode *Octree::find_closest_leaf(const Camera &cam)
{
	const Point3 &p = cam.origin();
	if (contains(p)) {
		OctreeNode *leaf = find_leaf_containing(p);
		assert(leaf != NULL);
		return leaf;
	}

	OctreeNode::Direction face_directions[3];
	int num_directions = 0;
	classify_point(p, face_directions, num_directions);
	OctreeNode *closest = NULL;
	real_t curr_dist = 0;
	real_t closest_dist = real::MAX;
	for (int i = 0; i < num_directions; i++) {
		OctreeNode::Direction d = face_directions[i];
		for (unsigned j = 0; j < _outside_nodes[d].size(); j++) {
			OctreeNode *node = _outside_nodes[d][j];
			if (!cam.sees(node->bounding_box())) {
				continue;
			}
			// @@@ FIXME: should be distance to face
			curr_dist = Point3::distance(
				p, node->centroid());
			if (curr_dist < closest_dist) {
				closest_dist = curr_dist;
				closest = node;
			}
		}
	}
	return closest;
}


void Octree::intersect(const Ray3 &ray, OctreeNode::hit_set &hits)
{
	assert(NULL != _root);
	_root->intersect(ray, hits);
}


int Octree::find_max_id() const
{
	assert(_root != NULL);
	return _root->max_id();
}


unsigned Octree::find_max_depth() const
{
	assert(_root != NULL);
	return _root->find_max_depth();
}


void Octree::init()
{
	compute_arrangement();
	init_occupancy_statistics();
	init_solidity();
}


// FIXME: we really need a scene graph with render context...

void Octree::render(const Camera &cam,
		    const ColorRgb &octree_color,
		    const ColorRgb &default_model_color,
		    bool show_model_enabled,
		    bool show_octree_enabled,
		    real_t lod,
		    bool points_as_cones) const
{
	assert(_root != NULL);
	unsigned num_points_rendered = 0;
	vector<const OctreeNode *> visible_leaves;
	collect_visible_leaves(cam, visible_leaves);
	for (unsigned i = 0; i < visible_leaves.size(); i++) {
		const OctreeNode *leaf = visible_leaves[i];
		assert(leaf != NULL);
		leaf->render(cam,
			     octree_color,
			     default_model_color,
			     show_model_enabled,
			     show_octree_enabled,
			     lod,
			     &num_points_rendered,
			     points_as_cones);
	}
	//fprintf(stderr, "num points rendered = %d\n", num_points_rendered);
}


bool Octree::has(OctreeNode *node) const
{
	assert(_root != NULL);
	return _root->has(node);
}


bool Octree::pick(const Ray3 &ray,
		  const char *model_name,
		  OctreeNode::hit_set &hits,
		  Point3 &picked_point,
		  OctreeNode *&picked_point_node)
{
	picked_point_node = NULL;
	bool picked = false;
	real_t min_hit_time = real::MAX;

	// find nodes hit by ray
	intersect(ray, hits);

	// for each node
	OctreeNode::hit_set::iterator hi;
	for (hi = hits.begin(); hi != hits.end(); hi++) {
		OctreeNode *node = hi->node();
		assert(node->is_leaf());

		if (node->num_vertices() == 0) {
			continue;
		}

		bool should_free_model = false;
		if (!node->has_model()) {
			node->read_data(model_name);
			should_free_model = true;
		}

		assert(node->has_model());

		// for each triangle in the node's model
		const IndexedTriangleSet &model = node->model();
		for (unsigned j = 0; j < model.num_triangles(); j++) {
			const IndexedTriangle &it = model.indexed_triangle(j);
			Polygon3 poly(model.vertex(it.A()),
				      model.vertex(it.B()),
				      model.vertex(it.C()));
			
			// intersect triangle with ray
			real_t hit_time;
			if (ray.intersects(poly, hit_time)) {
				picked = true;
				if (hit_time < min_hit_time) {
					min_hit_time = hit_time;
					picked_point_node = node;
				}
			}
		}

		if (should_free_model) {
			node->free_model();
		}
	}
	if (picked) {
		picked_point = ray.point(min_hit_time);
	}
	return picked;
}


bool Octree::pick_vertex(const Ray3 &ray,
			 real_t max_dist,
			 OctreeNode::hit_set &hits,
			 Point3 &picked_point)
{
	real_t min_dist = real::MAX;

	// find nodes hit by ray
	intersect(ray, hits);

	// for each node
	OctreeNode::hit_set::iterator hi;
	for (hi = hits.begin(); hi != hits.end(); hi++) {
		OctreeNode *node = hi->node();
		assert(node->is_leaf());

		if (node->num_vertices() == 0) {
			continue;
		}

		assert(node->has_model());

		// for each point in the node's model
		const IndexedTriangleSet &model = node->model();
		for (unsigned j = 0; j < model.num_vertices(); j++) {
			const Point3 &p = model.vertex(j);

			// check distance to ray
			real_t t = ray.t(p);
			if (real::is_positive(t)) {
				Point3 q = ray.point(t);
				real_t dist = Point3::distance(p, q);
				if (dist < min_dist) {
					picked_point = p;
					min_dist = dist;
				}
			}
		}
	}
	return min_dist < max_dist;
}


bool Octree::pick_vertex(const Ray3 &ray,
			 real_t max_dist,
			 OctreeNode::hit_set &hits,
			 Point3 &picked_point,
			 OctreeNodeCache &cache)
{
	real_t min_dist = real::MAX;

	// find nodes hit by ray
	intersect(ray, hits);

	// for each node
	OctreeNode::hit_set::iterator hi;
	for (hi = hits.begin(); hi != hits.end(); hi++) {
		OctreeNode *node = hi->node();
		assert(node->is_leaf());

		if (node->num_vertices() == 0) {
			continue;
		}

		if (!node->has_model()) {
			cache.fetch(node);
		}
		assert(node->has_model());

		// for each point in the node's model
		const IndexedTriangleSet &model = node->model();
		for (unsigned j = 0; j < model.num_vertices(); j++) {
			const Point3 &p = model.vertex(j);

			// check distance to ray
			real_t t = ray.t(p);
			if (real::is_positive(t)) {
				Point3 q = ray.point(t);
				real_t dist = Point3::distance(p, q);
				if (dist < min_dist) {
					picked_point = p;
					min_dist = dist;
				}
			}
		}
	}
	return min_dist < max_dist;
}


void Octree::print_leaves(const char *prefix) const
{
	assert(prefix != NULL);
	assert(_root != NULL);
	_root->print_leaves(prefix);
}


void Octree::scan_leaves_principal_components(const char *prefix)
{
	assert(_root != NULL);
	_root->scan_leaves_principal_components(prefix);
}


unsigned Octree::num_vertices() const
{
	assert(_root != NULL);
	unsigned n = 0;
	vector<OctreeNode *> leaves;
	_root->collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		n += leaf->num_vertices();
	}
	return n;
}


unsigned Octree::num_triangles() const
{
	assert(_root != NULL);
	unsigned n = 0;
	vector<OctreeNode *> leaves;
	_root->collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		OctreeNode *leaf = leaves[i];
		n += leaf->num_triangles();
	}
	return n;
}


unsigned Octree::find_max_num_triangles_per_leaf() const
{
	unsigned max_num_tri = 0;
	vector<const OctreeNode *> leaves;
	collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		const OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		if (leaf->num_triangles() > max_num_tri) {
			max_num_tri = leaf->num_triangles();
		}
	}
	return max_num_tri;
}


float Octree::find_avg_num_triangles_per_leaf() const
{
	float avg_num_tri = 0.0;
	vector<const OctreeNode *> leaves;
	collect_leaves(leaves);
	for (unsigned i = 0; i < leaves.size(); i++) {
		const OctreeNode *leaf = leaves[i];
		assert(leaf != NULL);
		avg_num_tri += ((float) leaf->num_triangles() /
				(float) leaves.size());
	}
	return avg_num_tri;
}


GTB_END_NAMESPACE
