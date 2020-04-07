
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


#ifndef __KDSPHERE_H
#define __KDSPHERE_H

/*
 * Author: Shachar Fleishman, shacharf@math.tau.ac.il
 *
 * Description:
 *   A KDTree that holds a set of spheres.
 *   The main purpose of this class is to cast a rays through spheres.
 *
 *   I assume the spheres are stored in some container with randon
 *   access indexing, i.e the [] operator is defined.
 *   In the tree only indices are stored.
 *
 *   Initialization:
 *
 *   Inserting spheres:
 *
 *   Casting a ray:
 *     xxx = InitializeRay(ray)
 *     while (xxx.CastRay(sphereList)
 *     {
 *         foreach sphere in sphereList - do something.
 *         if found intersection (or whatever condition) break;
 *     }
 *   
 * Template argument types:
 *    CONT - container type that holds the list of spheres.
 *
 * Changes:
 *
 * Who        When      What
 * ---------- --------- --------------------------------------------------
 * 
 */

#include "sphere.hpp"
#include "ray3.hpp"


template <class CONT>
class KDSphere
{
public:
    typdef unsigned index_type; // Type that is used to index into the
                                // user list of spheres.
    typedef std::vector spheres_list; // The output of this class is
                                // a list of spheres.
    class RayInfo;

    KDSphere(const CONT& cont);

    void insert(index_type idx);
    RayInfo InitializeRay(const ray3& ray);
    
    void MakeTree();

    /*
     * A class that holds information about a ray being cast.
     * More specifically, it holds a stack of cells that the ray
     * pierces. The stack grows and shrinks with time.
     *
     * The CastRay method is used by the user to retrieve a list
     * of spheres in the current-first cell that the ray pierces.
     */
    struct Node;
    class RayInfo
    {
    public:
        void insert(const Node* node);
        bool CastRay(spheres_list& result);

    private:
        typedef std::stack<const Node*> t_stack;
        t_stack _stack;
    };

private:
    const CONT& _spheres;

    //
    // Helper classes
    // 


    //
    // Node, parten to TreeNode and LeafNode, the structures
    // that hold the tree...
    //
    struct Node
    {
        Node() {}
        virtual ~Node() {}
        enum node_type {leaf, tree}; // Types of nodes

        virtual node_type Type() const=0;
    };

    struct TreeNode : public Node
    {
        ~TreeNode() { delete l; delete r; }
        node_type Type() const { return tree; }

        // Which axis is cut by this node
        // Numbering complies with VX, VY, VZ defined in algebra3.h
        enum t_axis {XAxis=0, YAxis, ZAxis}; 

        double cut_point;       // where along the axis is this node cut.
                                // Absolute coordinates.

        Node *l, *r;            // Left and right sons
        t_axis axis;
    }; // TreeNode

    struct LeafNode : public Node
    {
        LeafNode() {}
        node_type Type() const { return leaf; }
        
        void Insert(T v) { objects.push_back(v); }

        template<class IT>
        void Insert(IT first, IT last)
        {
            objects.insert(objects.end(), first, last);
        }

        //
        // Reserve space for additional K items, return previous # of elements
        //
        int Resize(int K)
        {
            int n = objects.size();
            objects.resize(n+K);
            return n;
        }

//        typedef tightvec<index_type> t_objects_list;
        typedef std::vector<index_type> t_objects_list;
        t_objects_list objects; // List of objects in the node

        typename t_objects_list::value_type& operator[](int n) {
                return objects[n];
        }
    }; // LeafNode

};

#endif // __KDSPHERE_H
