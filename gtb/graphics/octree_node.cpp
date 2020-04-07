
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
#include <gtb/graphics/octree_node.hpp>
#include <gtb/graphics/triangle3.hpp>
#include <gtb/graphics/camera.hpp>
#include <gtb/graphics/octree.hpp>
#include <gtb/io/io.hpp>
#include <gtb/error/error.hpp>
#include <gtb/util/util.hpp>
#include <map>
#include <set>
#include <string>
#endif // WIN32


#ifdef OUTLINE
#define inline
#include <gtb/graphics/octree_node.ipp>
#undef inline
#endif


GTB_BEGIN_NAMESPACE


using namespace std;


unsigned OctreeNode::_max_num_vertices = 500;
unsigned OctreeNode::_max_depth = 15;
unsigned OctreeNode::_global_ray_id = 0;


OctreeNode::OctreeNode(const Box3 &arg_box,
		       int arg_id,
		       Octant arg_octant,
		       unsigned arg_depth,
		       OctreeNode *arg_parent)
	: Box3(arg_box),
	  _id(arg_id),
	  _octant(arg_octant),
	  _depth(arg_depth),
	  _num_vertices(0),
	  _num_triangles(0),
	  _is_leaf(true),
	  _parent(arg_parent),
	  _model(NULL),
	  _ray_id(0)
{
	for(unsigned i = 0; i < NUM_CHILDREN; i++) {
		_children[i] = NULL;
	}
}


OctreeNode::~OctreeNode()
{
	for(unsigned i = 0; i < NUM_CHILDREN; i++) {
		if (_children[i] != NULL) {
			delete _children[i];
		}
	}
	if (_model != NULL) {
		delete _model;
	}
}


void OctreeNode::read_structure(FILE *fp)
{
	assert(fp != NULL);
	assert(sizeof(_octant) == 4);

	Box3::read(fp);
	read_int(&_id, fp);
    	read_unsigned((unsigned *) &_octant, fp);
    	read_unsigned(&_depth, fp);
    	read_unsigned(&_num_vertices, fp);
    	read_unsigned(&_num_triangles, fp);
    	read_bool(&_is_leaf, fp);

  	_parent = NULL;
  	for (unsigned i = 0; i < NUM_CHILDREN; i++) {
  		_children[i] = NULL;
  	}
  	_model = NULL;
	//fprintf(stderr, "read structure of node %d\n", _id);
}


void OctreeNode::write_structure(FILE *fp) const
{
	assert(fp != NULL);
	assert(sizeof(_octant) == 4);

	Box3::write(fp);
	write_int(_id, fp);
    	write_unsigned((unsigned) _octant, fp);
    	write_unsigned(_depth, fp);
    	write_unsigned(_num_vertices, fp);
    	write_unsigned(_num_triangles, fp);
    	write_bool(_is_leaf, fp);
	//fprintf(stderr, "wrote structure of node %d\n", _id);
}


void OctreeNode::read_data(const char *file_name_prefix)
{
	assert(is_leaf());
	assert(_num_vertices > 0);
	
	_model = Octree::read_data(file_name_prefix, _id);

//  	fprintf(stderr, "    read data of node %d\n", _id);
//  	fprintf(stderr, "    _num_vertices = %d\n", _num_vertices);
//  	fprintf(stderr, "    _model->num_vertices() = %d\n",
//  		_model->num_vertices());

	assert(_num_vertices == _model->num_vertices());
	assert(_num_triangles == _model->num_triangles());
}


void OctreeNode::write_data(const char *file_name_prefix) const
{
	assert(is_leaf());
	assert(_model != NULL);
	assert(_num_vertices > 0);

//  	fprintf(stderr, "    writing data of node %d\n", _id);
//  	fprintf(stderr, "    _num_vertices = %d\n", _num_vertices);
//  	fprintf(stderr, "    _model->num_vertices() = %d\n",
//  		_model->num_vertices());

	assert(_num_vertices == model().num_vertices());
	assert(_num_triangles == model().num_triangles());
	Octree::write_data(*_model, file_name_prefix, _id);
}


void OctreeNode::load_structure(FILE *fp, Octree &octree)
{
	read_structure(fp);
	if (!is_leaf()) {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			_children[i] = new OctreeNode(
				Box3(),	// box
				-1, // id
				ROOT, // octant
				0, // depth
				NULL); // parent
			assert(_children[i] != NULL);
			_children[i]->load_structure(fp, octree);
			_children[i]->_parent = this;
		}
	}
}


void OctreeNode::save_structure(FILE *fp) const
{
	write_structure(fp);
	if (!is_leaf()) {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(_children[i] != NULL);
			_children[i]->save_structure(fp);
		}
	}
}


void OctreeNode::load_data(const char *file_name_prefix)
{
	if (is_leaf()) {
		if (num_vertices() > 0) {
			read_data(file_name_prefix);
		}
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(_children[i] != NULL);
			_children[i]->load_data(file_name_prefix);
		}
	}
}


Box3 OctreeNode::child_box(unsigned i, unsigned j, unsigned k) const
{
	Point3 c = centroid();
	Point3 min_pt((i == 0) ? _min_pt[0] : c[0],
		      (j == 0) ? _min_pt[1] : c[1],
		      (k == 0) ? _min_pt[2] : c[2]);

	Point3 max_pt((i == 0) ? c[0] : _max_pt[0],
		      (j == 0) ? c[1] : _max_pt[1],
		      (k == 0) ? c[2] : _max_pt[2]);
	return Box3(min_pt, max_pt);
}


OctreeNode *OctreeNode::smallest_child_containing(const Box3 &box)
{
	if (!is_leaf()) {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			if (child(i)->classify_position(box) == INSIDE) {
				return child(i)->smallest_child_containing(
					box);
			}
		}
	}
	// We should really test for == INSIDE, but for now we have to
	// accept OVERLAP as well.  It may happen during an
	// out-of-core re-insertion.
	assert(classify_position(box) != OUTSIDE);
	return this;
}


void OctreeNode::collect_leaves(vector<const OctreeNode *> &leaves) const
{
	if (is_leaf()) {
		leaves.push_back(this);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->collect_leaves(leaves);
		}
	}
}


void OctreeNode::collect_leaves(vector<OctreeNode *> &leaves)
{
	if (is_leaf()) {
		leaves.push_back(this);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->collect_leaves(leaves);
		}
	}
}


void OctreeNode::collect_leaves_excluding_subtree(
	const OctreeNode *subtree_to_exclude,
	std::vector<OctreeNode *> &leaves)
{
	if (this == subtree_to_exclude) {
		return;
	} else if (is_leaf()) {
		leaves.push_back(this);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(_children[i] != NULL);
			_children[i]->collect_leaves_excluding_subtree(
				subtree_to_exclude, leaves);
		}
	}
}


void OctreeNode::collect_visible_leaves(
	const Camera &cam,
	vector<OctreeNode *> &nodes)
{
	// not visible
	if (!cam.sees(bounding_box())) {
		//fprintf(stderr, "node %d ignored\n", id());
		return;
	}

	// totally visible
	if (cam.contains(bounding_box())) {
		collect_leaves(nodes);
		return;
	}

	// partially visible
	if (is_leaf()) {
		nodes.push_back(this);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->collect_visible_leaves(cam, nodes);
		}
	}
}


void OctreeNode::collect_visible_leaves(
	const Camera &cam,
	vector<const OctreeNode *> &nodes) const
{
	// not visible
	if (!cam.sees(bounding_box())) {
		//fprintf(stderr, "node %d ignored\n", id());
		return;
	}

	// totally visible
	if (cam.contains(bounding_box())) {
		collect_leaves(nodes);
		return;
	}

	// partially visible
	if (is_leaf()) {
		nodes.push_back(this);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->collect_visible_leaves(cam, nodes);
		}
	}
}


void OctreeNode::collect_nodes(vector<OctreeNode *> &nodes)
{
	nodes.push_back(this);
	for (unsigned i = 0; i < NUM_CHILDREN; i++) {
		if (child(i) != NULL) {
			child(i)->collect_nodes(nodes);
		}
	}
}


unsigned OctreeNode::num_leaves() const
{
	if (is_leaf()) {
		return 1;
	} else {
		unsigned n = 0;
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			n += child(i)->num_leaves();
		}
		return n;
	}
}


unsigned OctreeNode::num_nodes() const
{
	if (is_leaf()) {
		return 1;
	} else {
		unsigned n = 0;
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			n += child(i)->num_nodes();
		}
		return n + 1;
	}
}


void OctreeNode::create_children(Octree &octree)
{
	for (unsigned i = 0; i < 2; i++) {
		for (unsigned j = 0; j < 2; j++) {
			for (unsigned k = 0; k < 2; k++) {
				Box3 b = child_box(i, j, k);

//  				fprintf(stderr, "%f %f %f %f %f %f\n",
//  					b.x_min(), b.y_min(), b.z_min(),
//  					b.x_max(), b.y_max(), b.z_max());

				Octant o = octant(i, j, k);
				_children[o] = new OctreeNode(
					b, // box
					octree.make_id(), // id
					o, // octant
					_depth + 1, // depth
					this); // parent
				assert(_children[o] != NULL);
			}
		}
	}
	_is_leaf = false;
	for (unsigned i = 0; i < NUM_CHILDREN; i++) {
		assert(child(i) != NULL);
	}
}


void OctreeNode::insert_vertices(const IndexedTriangleSet &its,
				 Octree &octree,
				 bool verbose,
				 FILE *log_file)
{
	assert(_vertices.size() == 0);
	timer t;
	t.start();
	if (verbose) {
		fprintf(stderr, "inserting vertices...");
	}
	for (unsigned i = 0; i < its.num_vertices(); i++) {
		if (!insert_vertex(i, its.vertices(), octree)) {
			GTB_ERROR("insert vertex failed");
		}
		if (verbose && (((i + 1) % 10000) == 0)) {
			fprintf(stderr, "\rinserting vertices... %d/%d",
				i + 1, its.num_vertices());
		}
	}
	if (verbose) {
		fprintf(stderr, "\rinserting vertices... %d/%d\n",
			its.num_vertices(), its.num_vertices());
	}
	t.stop();
	if (log_file != NULL) {
		fprintf(log_file, "insert vertices: %.3fs\n",
			t.elapsed_seconds());
	}
}


bool OctreeNode::insert_vertex(unsigned index,
			       const vector<Point3> &vertices,
			       Octree &octree)
{
	const Point3 &v = vertices[index];

	if (!is_leaf()) {
		return _children[octant(v)]->insert_vertex(
			index, vertices, octree);
	}

	_vertices.push_back(index);
	if (_vertices.size() <= _max_num_vertices) {
		return true;
	}

	if (_depth >= _max_depth) {
		GTB_ERROR("maximum octree depth reached");
		//return false;
	}

	create_children(octree);

	// distribute points
	for (unsigned m = 0; m < _vertices.size(); m++) {
		Point3 v2 = vertices[_vertices[m]];
		if (!_children[octant(v2)]->insert_vertex(
			_vertices[m], vertices, octree)) {
			return false;
		}
	}
	_vertices.clear();
	assert(_vertices.size() == 0);
	assert(_triangles.size() == 0);
	return true;
}


void OctreeNode::insert_triangles(const IndexedTriangleSet &its,
				  vector <unsigned> &replications,
				  bool verbose,
				  FILE *log_file)
{
	assert(_triangles.size() == 0);
	timer t;
	t.start();
	if (verbose) {
		fprintf(stderr, "inserting triangles...");
	}
	for (unsigned i = 0; i < its.num_triangles(); i++) {
		if (!insert_triangle(i, its.vertices(), its.triangles(),
				     replications)) {
			GTB_ERROR("insert triangle failed");
		}
		if (verbose && (((i + 1) % 10000) == 0)) {
			fprintf(stderr, "\rinserting triangles... %d/%d",
				i + 1, its.num_triangles());
		}
	}
	if (verbose) {
		fprintf(stderr, "\rinserting triangles... %d/%d\n",
			its.num_triangles(), its.num_triangles());
	}
	t.stop();
	if (log_file != NULL) {
		fprintf(log_file, "insert triangles: %.3fs\n",
			t.elapsed_seconds());
	}
}


bool OctreeNode::insert_triangle(unsigned index,
				 const vector<Point3> &vertices,
				 const vector<IndexedTriangle> &triangles,
				 vector<unsigned> &replications)
{
	// find smallest child containing triangle
	IndexedTriangle it = triangles[index];
	Box3 b = Box3::bounding_box(vertices[it.A()],
				    vertices[it.B()],
				    vertices[it.C()]);
	// It's OK to overlap.  It may happen during an out-of-core
	// reinsertion.
	assert(classify_position(b) != OUTSIDE);

	// recurse
	replications.push_back(0);
	insert_triangle(index, b, replications);

	return true;
}


void OctreeNode::insert_triangle(unsigned index,
				 const Box3 &b,
				 vector<unsigned> &replications)
{
	if (is_leaf()) {
		_triangles.push_back(index);
		replications[replications.size() - 1]++;
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			int pos = child(i)->classify_position(b);
			if (pos != OUTSIDE) {
				child(i)->insert_triangle(index, b,
							  replications);
				// early termination
				if (pos == INSIDE) {
					break;
				}
			}
		}
	}
}


void OctreeNode::create_submodel(const IndexedTriangleSet &its)
{
	assert(NULL == _model);
	assert((_vertices.size() > 0) || (_triangles.size() > 0));
	_model = its.create_submodel(_vertices, _triangles);
	assert(NULL != _model);
	_num_vertices = _model->num_vertices();
	_num_triangles = _model->num_triangles();
}


void OctreeNode::add_submodel(const IndexedTriangleSet &its)
{
	assert(NULL != _model);
	assert((_vertices.size() > 0) || (_triangles.size() > 0));
	IndexedTriangleSet *submodel = its.create_submodel(_vertices,
							   _triangles);
	assert(NULL != submodel);
	_model->add(*submodel);
	_num_vertices = _model->num_vertices();
	_num_triangles = _model->num_triangles();
}


void OctreeNode::free_model()
{
	assert(NULL != _model);
	delete _model;
	_model = NULL;
	_vertices.clear();
	_triangles.clear();
}


OctreeNode *OctreeNode::find_leaf_containing(const Point3 &p)
{
	if (is_leaf()) {
		assert(contains(p));
		return this;
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			if (child(i)->contains(p)) {
				OctreeNode *cell =
				    child(i)->find_leaf_containing(p);
				if (cell != NULL) {
					return cell;
				}
			}
		}
	}
	return NULL;
}


OctreeNode *OctreeNode::find_node(int node_id)
{
	if (_id == node_id) {
		return this;
	}
	if (!is_leaf()) {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			if (child(i)->find_node(node_id) != NULL) {
				return child(i);
			}
		}
	}
	return NULL;
}


bool OctreeNode::adj(unsigned d, unsigned o)
{
	bool ret = false;

	switch (d) {

	case L:
		switch (o) {
		case LDB:
		case LDF:
		case LUB:
		case LUF:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	case R:
		switch (o) {
		case RDB:
		case RDF:
		case RUB:
		case RUF:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	case D:
		switch (o) {
		case LDB:
		case LDF:
		case RDB:
		case RDF:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	case U:
		switch (o) {
		case LUB:
		case LUF:
		case RUB:
		case RUF:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	case B:
		switch (o) {
		case LDB:
		case LUB:
		case RDB:
		case RUB:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	case F:
		switch (o) {
		case LDF:
		case LUF:
		case RDF:
		case RUF:
			ret = true;
			break;
		default:
			ret = false;
			break;
		}
		break;

	default:
		ret = false;
		break;
	}
	return ret;
}


// Based on Hanan Samet, Applications of Spatial Data Structures
unsigned OctreeNode::reflect(unsigned dir, unsigned o)
{
	unsigned ret = ROOT;

	switch (dir) {
	case L:
		switch (o) {
		case LDB:
			ret = RDB;
			break;
		case LDF:
			ret = RDF;
			break;
		case LUB:
			ret = RUB;
			break;
		case LUF:
			ret = RUF;
			break;
		case RDB:
			ret = LDB;
			break;
		case RDF:
			ret = LDF;
			break;
		case RUB:
			ret = LUB;
			break;
		case RUF:
			ret = LUF;
			break;
		}
		break;
	case R:
		switch (o) {
		case LDB:
			ret = RDB;
			break;
		case LDF:
			ret = RDF;
			break;
		case LUB:
			ret = RUB;
			break;
		case LUF:
			ret = RUF;
			break;
		case RDB:
			ret = LDB;
			break;
		case RDF:
			ret = LDF;
			break;
		case RUB:
			ret = LUB;
			break;
		case RUF:
			ret = LUF;
			break;
		}
		break;
	case D:
		switch (o) {
		case LDB:
			ret = LUB;
			break;
		case LDF:
			ret = LUF;
			break;
		case LUB:
			ret = LDB;
			break;
		case LUF:
			ret = LDF;
			break;
		case RDB:
			ret = RUB;
			break;
		case RDF:
			ret = RUF;
			break;
		case RUB:
			ret = RDB;
			break;
		case RUF:
			ret = RDF;
			break;
		}
		break;
	case U:
		switch (o) {
		case LDB:
			ret = LUB;
			break;
		case LDF:
			ret = LUF;
			break;
		case LUB:
			ret = LDB;
			break;
		case LUF:
			ret = LDF;
			break;
		case RDB:
			ret = RUB;
			break;
		case RDF:
			ret = RUF;
			break;
		case RUB:
			ret = RDB;
			break;
		case RUF:
			ret = RDF;
			break;
		}
		break;
	case B:
		switch (o) {
		case LDB:
			ret = LDF;
			break;
		case LDF:
			ret = LDB;
			break;
		case LUB:
			ret = LUF;
			break;
		case LUF:
			ret = LUB;
			break;
		case RDB:
			ret = RDF;
			break;
		case RDF:
			ret = RDB;
			break;
		case RUB:
			ret = RUF;
			break;
		case RUF:
			ret = RUB;
			break;
		}
		break;
	case F:
		switch (o) {
		case LDB:
			ret = LDF;
			break;
		case LDF:
			ret = LDB;
			break;
		case LUB:
			ret = LUF;
			break;
		case LUF:
			ret = LUB;
			break;
		case RDB:
			ret = RDF;
			break;
		case RDF:
			ret = RDB;
			break;
		case RUB:
			ret = RUF;
			break;
		case RUF:
			ret = RUB;
			break;
		}
		break;
	default:
		GTB_ERROR("OctreeNode::reflect should never be here\n");
		break;
	}
	return ret;
}


OctreeNode *OctreeNode::neighbor(Direction dir)
{
	vector<unsigned> path;
	OctreeNode *o = this;

	while (1) {
		path.push_back(o->_octant);
		if (o->_parent == NULL) {
			return NULL;
		} else {
			if (!adj(dir, o->_octant)) {
				o = o->_parent;
				break;
			} else {
				o = o->_parent;
			}
		}
	}
	assert(path.size() > 0);
	for (int i = path.size() - 1; i > -1; i--) {
		unsigned ro = reflect(dir, path[i]);
		if (o->_depth == _depth) {
			break;
		}
		if (o->_children[ro] != NULL) {
			o = o->_children[ro];
		} else {
			return o;
		}
	}
	return o;
}


void OctreeNode::init_occupancy_statistics(int &min_num_tri,
					   int &max_num_tri,
					   int &total_num_tri,
					   int &total_num_nodes) const
{
	if (is_leaf()) {
		int num_tri = _num_triangles;
		if (num_tri < min_num_tri) {
			min_num_tri = num_tri;
		}
		if (num_tri > max_num_tri) {
			max_num_tri = num_tri;
		}
		total_num_tri += num_tri;
		total_num_nodes++;
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->init_occupancy_statistics(
				min_num_tri,
				max_num_tri,
				total_num_tri,
				total_num_nodes);
		}
	}
}


void OctreeNode::init_solidity(unsigned max_num_triangles)
{
	if (is_leaf()) {
		real_t s = ((real_t) _num_triangles /
			    (real_t) max_num_triangles);
		_initial_solidity = s;
		// _initial_solidity = s * s * s;
		// _initial_solidity = 0.0;
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->init_solidity(max_num_triangles);
		}
	}
}


//  static int test_nodes[] = {
//  	0, 5, 868, 876, 882, 889, 4692, 5025, 5124
//  };


//  static bool is_test_node(const OctreeNode *node)
//  {
//  	REQUIRE(node != NULL);
//  	unsigned n = sizeof(test_nodes) / sizeof(test_nodes[0]);
//  	for (unsigned i = 0; i < n; i++) {
//  		if (node->id() == test_nodes[i]) {
//  			return true;
//  		}
//  	}
//  	return false;
//  }


// FIXME: we really need a scene graph with render context...
void OctreeNode::render(const Camera &cam,
			const ColorRgb &octree_color,
			const ColorRgb &default_model_color,
			bool show_model_enabled,
			bool show_octree_enabled,
			real_t lod,
			unsigned *num_points_rendered,
			bool points_as_cones) const
{
	assert(is_leaf());
	if (show_model_enabled && (num_vertices() > 0)) {
		assert(has_model());
		model().render(cam,
			       default_model_color,
			       lod,
			       num_points_rendered,
			       points_as_cones);
	}
	if (show_octree_enabled) {
		octree_color.load();
		outline();
	}
}


void OctreeNode::render_solid_box() const
{
	glBegin(GL_QUADS);

	// back
	glNormal3f(0.0, 0.0, -1.0);
	glVertex3f(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _min_pt[2]);

	// front
	glNormal3f(0.0, 0.0, 1.0);
	glVertex3f(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _max_pt[2]);

	// left
	glNormal3f(-1.0, 0.0, 0.0);
	glVertex3f(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _min_pt[2]);

	// right
	glNormal3f(1.0, 0.0, 0.0);
	glVertex3f(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _max_pt[2]);

	// top
	glNormal3f(0.0, 1.0, 0.0);
	glVertex3f(_min_pt[0], _max_pt[1], _min_pt[2]);
	glVertex3f(_min_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _max_pt[2]);
	glVertex3f(_max_pt[0], _max_pt[1], _min_pt[2]);

	// bottom
	glNormal3f(0.0, -1.0, 0.0);
	glVertex3f(_min_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _min_pt[2]);
	glVertex3f(_max_pt[0], _min_pt[1], _max_pt[2]);
	glVertex3f(_min_pt[0], _min_pt[1], _max_pt[2]);

	glEnd();
}


void OctreeNode::render_geometry() const
{
	model().render_geometry();
}


Vector3 OctreeNode::normal(Direction d)
{
	Vector3 v(Vector3::VECTOR3_ZERO);

	switch(d) {
	case L:
		v = Vector3(1, 0, 0);
		break;
	case R:
		v = Vector3(-1, 0, 0);
		break;
	case U:
		v = Vector3(0, -1, 0);
		break;
	case D:
		v = Vector3(0, 1, 0);
		break;
	case F:
		v = Vector3(0, 0, -1);
		break;
	case B:
		v = Vector3(0, 0, 1);
		break;
	default:
		fprintf(stderr, "%s:%d: invalid direction\n",
			__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return v;
}


int OctreeNode::max_id() const
{
	if (is_leaf()) {
		return id();
	} else {
		int maxid = -1;
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			int tmpid = child(i)->max_id();
			if (tmpid > maxid) {
				maxid = tmpid;
			}
		}
		assert(maxid >= 0);
		return maxid;
	}
}


unsigned OctreeNode::find_max_depth() const
{
	if (is_leaf()) {
		return depth();
	} else {
		unsigned maxdepth = 0;
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			unsigned tmpdepth = child(i)->find_max_depth();
			if (tmpdepth > maxdepth) {
				maxdepth = tmpdepth;
			}
		}
		return maxdepth;
	}
}


bool OctreeNode::has(OctreeNode *node) const
{
	if (is_leaf()) {
		if (node == this) {
			return true;
		}
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			if (child(i)->has(node)) {
				return true;
			}
		}
	}
	return false;
}


void OctreeNode::intersect(const Ray3 &ray, hit_set &hits)
{
	hits.clear();
	increment_global_ray_id();
	rec_intersect(ray, hits);
}


void OctreeNode::rec_intersect(const Ray3 &ray, hit_set &hits)
{
	// check if already visited
	if (_ray_id == _global_ray_id) {
		return;
	}
	_ray_id = _global_ray_id;

	// check bounding box
	real_t t1, t2;
	if (!ray.intersects(bounding_box(), t1, t2)) {
		return;
	}

	if (is_leaf()) {
		// add itself
		unsigned n = hits.size();
		hits.insert(hit(this, t1, t2));
		if (hits.size() != (n + 1)) {
			GTB_ERROR("hits.insert failed: check for duplicates");
		}
		// recurse on neighbors
		for (unsigned i = 0; i < _neighbors.size(); i++) {
			assert(_neighbors[i] != this);
			_neighbors[i]->rec_intersect(ray, hits);
		}
	} else {
		// recurse on children
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->rec_intersect(ray, hits);
		}
	}
}


void OctreeNode::print_leaves(const char *prefix) const
{
	assert(prefix != NULL);
	if (is_leaf()) {
		assert(_model != NULL);
		if (_model->num_vertices() < MIN_PCA_POINTS) {
			return;
		}
		char id_str[20];
		sprintf(id_str, "%d", _id);
		string file_name = (string(prefix) +
				    string("-") +
				    string(id_str) +
				    string(".leaf"));
		//fprintf(stderr, "file name=%s\n", file_name.c_str());
		FILE *fp = fopen(file_name.c_str(), "wb");
		if (fp == NULL) {
			perror(file_name.c_str());
			exit(EXIT_FAILURE);
		}
		fprintf(fp, "%d\n", _model->num_vertices());
		for (unsigned i = 0; i < _model->num_vertices(); i++) {
			const Point3 &v = _model->vertex(i);
			fprintf(fp, "%f %f %f\n", v.x(), v.y(), v.z());
		}
		fclose(fp);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->print_leaves(prefix);
		}
	}
}


void OctreeNode::scan_leaves_principal_components(const char *prefix)
{
	assert(prefix != NULL);
	if (is_leaf()) {
		assert(_model != NULL);
		if (_model->num_vertices() < MIN_PCA_POINTS) {
			return;
		}
		char id_str[20];
		sprintf(id_str, "%d", _id);
		string file_name = (string(prefix) +
				    string("-") +
				    string(id_str) +
				    string(".leaf.pca"));

		FILE *fp = fopen(file_name.c_str(), "rb");
		if (fp == NULL) {
			perror(file_name.c_str());
			exit(EXIT_FAILURE);
		}

		float o1, o2, o3;
		float x1, x2, x3;
		float y1, y2, y3;
		float z1, z2, z3;
		float l1, l2, l3;

		// read origin
		fscanf(fp, "%f %f %f", &o1, &o2, &o3);

		// read axes: by convention, y is the first principal component
		fscanf(fp, "%f %f %f", &y1, &y2, &y3);
		fscanf(fp, "%f %f %f", &x1, &x2, &x3);
		fscanf(fp, "%f %f %f", &z1, &z2, &z3);
		fscanf(fp, "%f %f %f", &l1, &l2, &l3);

		Point3 o(o1, o2, o3);
		Vector3 x = (real_t)l2 / l3 * Vector3(x1, x2, x3);
		Vector3 y = (real_t)l1 / l3 * Vector3(y1, y2, y3);
		Vector3 z = Vector3(z1, z2, z3);
        if (y.dot(o - Point3::ZERO) > (real_t)0.0) {
			y = -y;
		}
		_cs.reset(o, x, y, z);

		fclose(fp);
	} else {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->scan_leaves_principal_components(prefix);
		}
	}
}


void OctreeNode::increment_depth()
{
	if (_depth >= _max_depth) {
		GTB_ERROR("maximum octree depth reached");
	}
	_depth++;
	if (!is_leaf()) {
		for (unsigned i = 0; i < NUM_CHILDREN; i++) {
			assert(child(i) != NULL);
			child(i)->increment_depth();
		}
	}
}


void OctreeNode::adopt(OctreeNode *kid, Octree &octree)
{
	assert(kid != NULL);
	assert(kid->_depth == 0);
	assert(kid->_octant == ROOT);
	assert(_depth == 0);
	assert(_octant == ROOT);

	// get octant of adopted kid
	assert(contains(kid->centroid()));
	OctreeNode::Octant o1 = octant(kid->centroid());
	kid->_octant = o1;

	// adopt kid
	assert(_children[o1] == NULL);
	_children[o1] = kid;
	kid->_parent = this;
	_is_leaf = false;

	// update kid's depth
	kid->increment_depth();

	// create empty siblings for the adopted kid
	for (unsigned i = 0; i < 2; i++) {
		for (unsigned j = 0; j < 2; j++) {
			for (unsigned k = 0; k < 2; k++) {
				Octant o2 = octant(i, j, k);
				if (o2 == o1) {
					continue;
				}
				assert(_children[o2] == NULL);
				Box3 b = child_box(i, j, k);
				_children[o2] = new OctreeNode(
					b, // box
					octree.make_id(), // id
					o2, // octant
					_depth + 1, // depth
					this); // parent
			}
		}
	}

	// sanity checks
	assert(_children[o1] == kid);
	for (int i = 0; i < NUM_CHILDREN; i++) {
		assert(_children[i] != NULL);
		if (i == o1) {
			assert(_children[i] == kid);
		} else {
			assert(_children[i] != kid);
		}
		assert(_children[i]->_parent == this);
	}
	assert(kid->_depth == 1);
	assert(kid->_octant != ROOT);
	assert(_depth == 0);
	assert(_octant == ROOT);
}


GTB_END_NAMESPACE
