
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


#ifndef __CONEKDTREE_H
#define __CONEKDTREE_H


class Box4 {
    public:
    Box4() {
    }
    Box4(Box3 _b, real_type rmin, real_type rmax) {
	b = _b;
	r[0] = rmin;
	r[1] = rmax;
    }

    Box3 b;
    real_type r[2];
};

inline Box4 make_union(const Box4 &lhs, const Box4 &rhs) {
    Box4 ret;
    ret.b = Box3::make_union(lhs.b, rhs.b);
    ret.r[0] = std::min(lhs.r[0], rhs.r[0]);
    ret.r[1] = std::max(lhs.r[1], rhs.r[1]);
    return ret;
}

class Cone {
    public:
    Point3 p;
    real_type r;
    real_type angle_dot;
};

inline bool intersect(const Box4 &b, const Cone &c) {

    if (c.r > b.r[1])
	return false;

    if (b.b.x_min() == c.p[0] &&b.b.x_max() == c.p[0] &&
	b.b.y_min() == c.p[1] &&b.b.y_max() == c.p[1] &&
	b.b.z_min() == c.p[2] &&b.b.z_max() == c.p[2] &&
	b.r[0] == c.r && b.r[1] == c.r)
	return false;
	

    if (b.b.contains(c.p)) {
	return true;
    } else {

	Point3 closest = c.p;
	if (closest[0] < b.b.x_min())   closest[0] = b.b.x_min();
	if (closest[0] > b.b.x_max())   closest[0] = b.b.x_max();
	if (closest[1] < b.b.y_min())   closest[1] = b.b.y_min();
	if (closest[1] > b.b.y_max())   closest[1] = b.b.y_max();
	if (closest[2] < b.b.z_min())   closest[2] = b.b.z_min();
	if (closest[2] > b.b.z_max())   closest[2] = b.b.z_max();
	real_type r = b.r[1];

	Vector3 diff = closest - c.p;
	real_type diff_r = r - c.r;

	real_type len = sqrt(diff.dot(diff) + diff_r*diff_r);
	diff /= len;
	diff_r /= len;

	real_type thisdot = diff_r;

	return (thisdot > c.angle_dot);
    }
}




// DISTBOXCLASS	must provide two typenames and two functions:
// typename DISTBOXCLASS::Box3
// typename DISTBOXCLASS::Point3
// DISTBOXCLASS::bounding_box(const OBJECT &) const;
// DISTBOXCLASS::distance(const OBJECT &, const typename DISTBOXCLASS::Point3 &) const;

template <typename OBJECT, typename DISTBOXCLASS>
class ConeBoxKDTree {

    public:

#define CONEBOXKDTREE_MAXLEAFSIZE	10

    ConeBoxKDTree() {
	children[0] = children[1] = NULL;
	intersected = false;
    }
    ~ConeBoxKDTree() {
	if (children[0]) delete children[0];
	if (children[1]) delete children[1];
    }


    ConeBoxKDTree(const std::vector<OBJECT> &iobjects, const DISTBOXCLASS &boxclass, int axis=0) {
	intersected = false;

	children[0] = children[1] = NULL;

	ReBuild(iobjects, boxclass, axis);
    }

    void ReBuild(const std::vector<OBJECT> &iobjects, const DISTBOXCLASS &boxclass, int axis=0) {

	if (children[0]) delete children[0];
	if (children[1]) delete children[1];
	children[0] = children[1] = NULL;

	objects = iobjects;

	// make the bounding box
        bbox = boxclass.bounding_box(objects[0]);
	for (unsigned i=1; i<objects.size(); i++) {

	    bbox = make_union(bbox, boxclass.bounding_box(objects[i]));
	}

	Split(boxclass, axis);

    }


    // get all the objects who's bounding boxes intersect the given cone
    bool MarkIntersected(const DISTBOXCLASS &boxclass, const Cone &c, vector<int> &marked) {

	if (intersected)
	    return true;

	// check our bounding box
	if (!intersect(bbox, c)) {
	    return false;
	}

	// check any leaf objects
	bool all_children_intersected = true;
	for (unsigned i=0; i<objects.size(); i++) {
	    Box4 obox = boxclass.bounding_box(objects[i]);

	    if (intersect(obox, c)) {
		marked[objects[i]] = 1;
	    } else {
		all_children_intersected &= (marked[objects[i]]==1);
	    }
	}


	// try going into the children
	bool all_leafs_intersected = true;
	if (children[0])	all_leafs_intersected &= children[0]->MarkIntersected(boxclass, c, marked);
	if (children[1])	all_leafs_intersected &= children[1]->MarkIntersected(boxclass, c, marked);


	if (all_children_intersected && all_leafs_intersected)
	    intersected = true;

	return intersected;
    }



    void Split(const DISTBOXCLASS &boxclass, int axis=-1) {


	// check if we should stop splitting
	if (objects.size() >= CONEBOXKDTREE_MAXLEAFSIZE) {

	    // if we're not told what axis to use, use the biggest
	    if (axis == -1) {
		real_type xl, yl, zl, rl;
		xl = bbox.b.x_length();
		yl = bbox.b.y_length();
		zl = bbox.b.z_length();
		rl = bbox.r[1]-bbox.r[0];

		if (xl>yl && xl>zl && xl>rl) axis = 0;
		else if (yl>xl && yl>zl && yl>rl) axis = 1;
		else if (zl>xl && zl>yl && zl>rl) axis = 2;
		else axis = 3;
	    }


	    // split the list by the axis
	    std::vector<OBJECT> cobjects[2];

	    if (objects.size() < 500) {
		std::vector< std::pair<real_type,OBJECT> > sorter(objects.size());

		for (unsigned i=0; i<objects.size(); i++) {
		    Box4 obox = boxclass.bounding_box(objects[i]);
		    real_type sortkey;
		    switch (axis) {
		    case 0:
			sortkey = obox.b.x_max()+obox.b.x_min(); break;
		    case 1:
			sortkey = obox.b.y_max()+obox.b.y_min(); break;
		    case 2:
			sortkey = obox.b.z_max()+obox.b.z_min(); break;
		    case 3:
			sortkey = obox.r[1]+obox.r[0]; break;
		    }
		    sortkey /= 2;

		    sorter[i] = std::pair<real_type,OBJECT>(sortkey, objects[i]);
		}

		std::sort(sorter.begin(), sorter.end());


		unsigned i;
		for (i=0; i<sorter.size()/2; i++) {
		    cobjects[0].push_back(sorter[i].second);
		}
		for ( ; i<sorter.size(); i++) {
		    cobjects[1].push_back(sorter[i].second);
		}

	    } else {
		for (unsigned i=0; i<objects.size(); i++) {
		    Box4 obox = boxclass.bounding_box(objects[i]);

		    real_type sortkey;
		    switch (axis) {
		    case 0:
			sortkey = obox.b.x_max()+obox.b.x_min(); break;
		    case 1:
			sortkey = obox.b.y_max()+obox.b.y_min(); break;
		    case 2:
			sortkey = obox.b.z_max()+obox.b.z_min(); break;
		    case 3:
			sortkey = obox.r[1]+obox.r[0]; break;
		    }
		    sortkey /= 2;


		    if (axis < 3) {
			if (sortkey < bbox.b.centroid()[axis]) {
			    cobjects[0].push_back(objects[i]);
			} else {
			    cobjects[1].push_back(objects[i]);
			}
		    } else {
			if (sortkey < (bbox.r[0]+bbox.r[1])/2) {
			    cobjects[0].push_back(objects[i]);
			} else {
			    cobjects[1].push_back(objects[i]);
			}
		    }
		}
	    }


	    if ((cobjects[0].size() != 0) && (cobjects[1].size() != 0)) {

		// actually have to split

		objects.clear();
		compact(objects);

		assert(!children[0] && !children[1]);

		children[0] = new ConeBoxKDTree(cobjects[0], boxclass, (axis+1)%4);
		children[1] = new ConeBoxKDTree(cobjects[1], boxclass, (axis+1)%4);

	    }
	}
    }




    // leaf node
    std::vector<OBJECT> objects;

    // internal node
    ConeBoxKDTree* children[2];

    Box4  bbox;
    bool intersected;

};








#endif // __CONEKDTREE_H
