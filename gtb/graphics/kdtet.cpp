
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
#include "kdtet.h"

GTB_BEGIN_NAMESPACE

KDTET::KDTET(int max_objects_in_node) :
    _max_objects_in_node(max_objects_in_node),
    _tree(0)
{
    _tree = new LeafNode;
}

//
// BUGBUG: support for empty tree only
//
void KDTET::Insert(p_object pobj)
{
    assert(_tree->type == Node::leaf);
    LeafNode* rootleaf = static_cast<LeafNode*>(_tree);

    rootleaf->Insert(pboj);
}


void KDTET::MakeTree()
{
    Node* root = _tree.release();
    _tree = root->MakeTree(_max_objects_in_node);
}

p_object KDTET::FindTet(const Point3& p)
{
    Node* root = _tree;

    while (root->type == Node::tree)
    {
        TreeNode* treenode = root;
        if (p[treenode->axis] < treenode->cut_point) root = treenode->l;
        else root = treenode->r;
        assert(root != 0);
    }

    assert(root->type == Node::leaf);
    LeafNode* leaf = root;

    p_object result;
    if (leaf->Contains(p, result)) return result;
    else return 0;
}

Node* KDTET::TreeNode::MakeTree(int max_objects_in_leaf)
{
    l = l->MakeTree(max_objects_in_leaf);
    r = r->MakeTree(max_objects_in_leaf);
    return this;
}

Node* KDTET::LeafNode::MakeTree(int max_objects_in_leaf)
{
    if (_objects.size() <= max_objects_in_leaf) return this;

    //
    // Compute a split point and axis with the following heuristic:
    // Choose three points at random
    // compute their bounding box
    // Take the split axis to be the maximum of the width,height,depth
    //   of the bounding box
    // Take the position of the split at the middle of the bounding box
    //
    TreeNode::t_axis axis;
    double cut_point;
    {
        unsigned idx1 = rand32() % n_objects;
        unsigned idx2 = rand32() % n_objects;
        unsigned idx3 = rand32() % n_objects;
        const Point3& p1 = _objects[idx1]->centroid();
        const Point3& p2 = _objects[idx2]->centroid();
        const Point3& p3 = _objects[idx3]->centroid();

        Point3 min_point(
            min3(p1[0],p2[0],p3[0]),
            min3(p1[1],p2[1],p3[1]),
            min3(p1[2],p2[2],p3[2]));
        Point3 max_point(
            max3(p1[0],p2[0],p3[0]),
            max3(p1[1],p2[1],p3[1]),
            max3(p1[2],p2[2],p3[2]));

        Vector3 box_dimentions = max_point - min_point;

        int k;
        max3(box_dimentions[0], box_dimentions[1], box_dimentions[2], k);
        //                double an_epsilon = box_dimentions[k]/100;
        axis = static_cast<TreeNode::t_axis>(k);
        cut_point = min_point[axis] + box_dimentions[axis]*0.51;
    }

    //
    // Make the split
    //
    LeafNode* rightnode = new LeafNode;

    typename t_objects_list::iterator last_left = _objects.end();
    typename t_objects_list::iterator f = _objects.begin();
    while (f < last_left)
    {
        if (f->centroid()[axis] <= cut_point)
        {
            ++f;
        }
        else
        {
            --last_left;
            rightnode->_objects.push_back(*f);
            *f = *last_left;
        }
    }

    assert(last_left != _objects.end());

    _objects.erase(last_left, cell->objects.end());

    TreeNode* subroot = new TreeNode;
    subroot->cut_point = cut_point;
    subroot->axis = axis;
    subroot->l = MakeTree(max_objects_in_leaf);
    subroot->r = rightnode->MakeTree(max_objects_in_leaf);
    return subroot;
}

bool KDTET::LeafNode::Contains(const Point3& p, p_object& tet)
{
    typename t_objects_list::iterator f = _objects.begin();
    typename t_objects_list::iterator l = _objects.end();
    for (; f != l; ++f)
    {
        if (f->contains(p))
        {
            tet = *f;
            return true;
        }
    }
    tet = 0;
    return false;
}


GTB_END_NAMESPACE
