
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


#ifndef __KDTREE_H
#define __KDTREE_H

#include <gtb/stlext.h>
#include <gtb/graphics/box3.hpp>

/*
 * Author: Shachar Fleishman, shacharf@math.tau.ac.il
 * Description:
 *
 * A (hopefully) stand-alone 3D kd-tree class
 * Assumptions:
 *   1. uses the algebra3 vector objects
 *   2. Objects are points.
 *   3. There is a function GetPoint(object) that returns a reference
 *      to a tPoint3<REAL>. (Assumed to be cheep).
 *      replace by template argument
 *
 * Supported features:
 *  - Insert objects. (does not update the tree)
 *    NOTE: the implementation is memory efficient, inserting
 *          one object at a time is slow, it is better
 *          to use the insert(f,l) function
 *
 *  - MakeTree - Computes the tree.
 *    Insert, MakeTree can be interchanged, i.e.
 *    Insert, MakeTree loop.
 *  - user-defined minimal # of objects/cell
 *  - Tranverse objects in the neighborhood of a point to a given radius / number.
 *  - Traversal methods:
 *    1. by giving a function that will be called with every object
 *       that is in the neighborhood
 *    2. init / getnext mechanism
 *
 * Future:
 *  - Ray trace through tree.
 *  - Clean up traversal and insert
 *  - Make fast traversal non-ordered for points in radius
 *  - Can I make quick epsilon accuracy traversal?
 *
 *  - traverse search with distance^2!!!
 *  - Remove the owner = false in traversel, check scalable???
 *  - remove the points color from render!!!
 */

// SHACHAR: I think there is a bug in the find / traverse
// routines, it should be 
// ..
//     if (d >= 0) // Mark2
//     {
// ...
//
// and not
//
// ..
//     if (d > 0) // Mark2
//     {
// ...
//
//

//#include <assert.h>
#include <vector>
#include <queue>
#include "vector3.hpp"
#include "point3.hpp"
//#include <tightvec.h>
//#include <tools.h>
#include <gtb/memory/ptrs.h>
#include <gtb/memory/hirestimer.h>

GTB_BEGIN_NAMESPACE

#define __KDT_PROFILING 0

#define FASTPQ2 0

//
// Profiling
//
extern CHighResTimer _timer_push;
extern CHighResTimer _timer_pop;
extern CHighResTimer _timer_insert;
extern CHighResTimer _timer_makeheap;
extern CHighResTimer _timer_ordered_init_pq;
extern CHighResTimer _timer_ordered_init_pq2;
extern CHighResTimer _timer_ordered_getnext;

#define PROFILE_UNORDERED_QUERY 0
#define DEBUG_TRAVERSE 0

struct kdtree_exception_no_more_values {};
extern kdtree_exception_no_more_values _exception_no_more_values;

#define KDT_INFINIT_NEIGHBORS 0x7fffffff

/*-------------------------- My priority queue --------------*/
/*
 * A priority queue that allows the user to insert many
 * objects before generating the pq.
 * Speed up...
 */
template<
    class T,
    class Cont = std::vector<T>,
    class Pred = std::less<typename Cont::value_type> >
class fast_pq
{
public:
    typedef typename Cont::value_type      value_type;
    typedef typename Cont::size_type       size_type;
    typedef          Cont                  container_type;

    typedef typename Cont::reference       reference;
    typedef typename Cont::const_reference const_reference;
protected:
    Cont c;
    Pred _comp;
public:
    fast_pq() : c() {c.reserve(100);}
    explicit fast_pq(const Pred& __x) :  c(), _comp(__x) {c.reserve(100);}
    fast_pq(const Pred& __x, const Cont& __s) 
        : c(__s), _comp(__x)
    { 
        make_heap(c.begin(), c.end(), _comp); 
    }

    //
    // standard priority queue routines
    //
    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    void clear() { c.clear(); }
    const_reference top() const { return c.front(); }

    void push(const value_type& __x) 
    {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_push);
#endif // __KDT_PROFILING==1
            c.push_back(__x); 
            std::push_heap(c.begin(), c.end(), _comp);
    }
    void pop() {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_pop);
#endif // __KDT_PROFILING==1
            std::pop_heap(c.begin(), c.end(), _comp);
            c.pop_back();
    }

    // My extentions

    template <class _InputIterator>
    void insert( _InputIterator __first, _InputIterator __last)
    {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_insert);
#endif // __KDT_PROFILING==1
        c.insert(c.end(), __first, __last);
        std::make_heap(c.begin(), c.end(), _comp);
    }  

    //
    // insert an element without remaking the heap
    // it is the users responsiblity to remake the heap
    //
    // BUG: should have added an assert on top(), pop()
    // in case the user did not call make heap
    //
    void push_only(const T& x)
    {
        c.push_back(x);
    }

    //
    //
    //
    void remake_heap()
    {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_makeheap);
#endif // __KDT_PROFILING==1
        std::make_heap(c.begin(), c.end(), _comp);
    }

    //
    // DEBUG:
    //    print
    //
    void print()
    {
        typename Cont::const_iterator f = c.begin();
        typename Cont::const_iterator l = c.end();
        for (; f != l; ++f)
        {
            printf("   "); f->print_type(); printf(" %g\n", f->GetDist());
        }
    }

}; // fast_pq

/*---------------------- [ My priority queue ] --------------*/

/*-------------------------- My priority queue 2 --------------*/
/*
 * tree-based pq.
 */
template<class T>
class fast_pq2
{
    typedef T& reference;
    typedef const T& const_reference;
    typedef unsigned size_type;
protected:
    std::multiset<T> c;
public:
    fast_pq2() : c() {}

    //
    // standard priority queue routines
    //
    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    const_reference top() const { return *c.begin(); }

    void push(const T& _x) 
    {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_push);
#endif // __KDT_PROFILING==1
        c.insert(_x);
    }
    void pop() {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_pop);
#endif // __KDT_PROFILING==1
        c.erase(c.begin());
    }


    //
    // insert an element without remaking the heap
    // it is the users responsiblity to remake the heap
    //
    // BUG: should have added an assert on top(), pop()
    // in case the user did not call make heap
    //
    void push_only(const T& x)
    {
        push(x);
    }

    //
    //
    //
    void remake_heap()
    {
    }

    //
    // DEBUG:
    //    print
    //
    void print()
    {
        typename std::multiset<T>::const_iterator f = c.begin();
        typename std::multiset<T>::const_iterator l = c.end();
        for (; f != l; ++f)
        {
            printf("   "); f->print_type(); printf(" %g\n", f->GetDist());
        }
    }
}; // fast_pq2

/*---------------------- [ My priority queue2 ] --------------*/


/*------------   Default GetPoint function object   ---------*/
template <class REAL>
struct t_default_get_point
{
    const tPoint3<REAL>& operator() (const tPoint3<REAL>& p) const { return p; }
};

template <class REAL>
t_default_get_point<REAL> gen_default_get_point()
{
    t_default_get_point<REAL> obj;
    return obj;
}

/*------------ [ Default GetPoint function object ] ---------*/

/*---------------------- Helper functions --------------*/
template <class REAL>
void SplitBBox(
    const tBox3<REAL>& bbox,      // Input bounding box
    int axis,
    REAL position,        // axis, position define where to split
    tBox3<REAL>& lbbox, 
    tBox3<REAL>& rbbox            // Result is store in [lr]bbox
    );
/*---------------------- [ Helper functions ] --------------*/


//
// A function object that compares two elements (vectors) along
// the desired axis
//
// Defined outside of KDTree because there is no trinary_function in stl.
//
// template arguments:
//    T - the type of object that is stored in the tree
//    GETPOINT - a function T ==> tPoint3<REAL>
//
template<class T, class GETPOINT>
struct KDTObjLess : public std::binary_function<T,T,bool>
{
    int axis;
    GETPOINT _getpoint;
    
    KDTObjLess(int AXIS, GETPOINT& getpoint) : axis(AXIS), _getpoint(getpoint) {}

    bool operator() (const T& lhs, const T&rhs)
    {
        return _getpoint(lhs)[axis] < _getpoint(rhs)[axis];
    }
};

template<class T, class GETPOINT>
inline KDTObjLess<T, GETPOINT> GetKDTObjLess(int axis, GETPOINT getpoint)
{
    return KDTObjLess<T, GETPOINT>(axis, getpoint);
}

//
// template arguments:
//    T - the type of object that is stored in the tree
//    REAL - type of point, i.e. tPoint3<REAL><REAL> or...
//    GETPOINT - a function T ==> tPoint3<REAL>
//
template<class T, class REAL, class GETPOINT = t_default_get_point<REAL> >
class KDTree
{
public:
    typedef unsigned int size_type;
    typedef std::vector<T> result_sequence;

    KDTree(size_type MAX_IN_CELL, const tBox3<REAL>& BBOX, const GETPOINT& getpoint = gen_default_get_point<REAL>());
    ~KDTree() { delete root; }

    const tBox3<REAL>& GetBBox() const { return bbox; }
    size_type GetMaxInCell() const { return max_in_cell; }

     void lprintf(const char *fmt, ...) const
     {
#if 0
	     va_list args;
	     va_start(args, fmt);
	     vfprintf(stdout, fmt, args);
	     va_end(args);
#endif
     }

    //
    // Insert an object.
    // if split is true, split the tree if necessary
    //
    void Insert(T object, bool split=false);
    void Remove(T object);	// JS


    //
    // Find the closet point
    //
    struct SaveSingle
    {
        mutable T p;
        void operator() (const T& v) const { p = v; }
    };

    T FindMin(const tPoint3<REAL>& p) const;

    //
    // Debug print
    // Requires function print(T);
    //
    void Print() const { Print(root, bbox, 0); }
    void PrintStatistics() const;

// Documentation: implementation is bellow...
//    template <class IT>
//    void Insert(IT first, IT last)

    void MakeTree();
    

    // Info only, defined bellow
    //template <class FUNC> void Traverse(
    //    const tPoint3<REAL>& point, 
    //    size_type K,
    //    REAL radius, // Maximal radius
    //    FUNC op);

    //
    // Traverse:
    //   The main method for extracting a set of ponits
    //   traverse calls op() for each x_i that is near the input point
    //
    template <class FUNC> void Traverse(const tPoint3<REAL>& point, size_type K, const FUNC& op)
    {
        Traverse(point, K, 1e38, op);
    }

    template <class FUNC> void Traverse(const tPoint3<REAL>& point, REAL radius, const FUNC& op)
    {
        Traverse(point, KDT_INFINIT_NEIGHBORS, radius, op);
    }

    template <class FUNC, class PRED> void Traverse_if(const tPoint3<REAL>& point, size_type K, const FUNC& op, const PRED& pred)
    {
        Traverse_if(point, K, 1e38, op, pred);
    }

    template <class FUNC, class PRED> void Traverse_if(const tPoint3<REAL>& point, REAL radius, const FUNC& op, const PRED& pred)
    {
        Traverse_if(point, KDT_INFINIT_NEIGHBORS, radius, op, pred);
    }

    //
    // A common use for traverse is to extract a set of points
    //
    template <class ITERATOR>
    void Extract(const tPoint3<REAL>& point, size_type K, REAL radius, ITERATOR it) const
    {
        Traverse(point, K, radius, AssignOp<T>(it));
    }

    template <class ITERATOR>
    void UnorderedExtract(const tPoint3<REAL>& point, REAL radius, ITERATOR it)
    {
        UnorderedTraverse(point, radius, AssignOp<T>(it));
    }

    template <class ITERATOR, class PRED>
    void Extract_if(const tPoint3<REAL>& point, size_type K, REAL radius, ITERATOR it, const PRED& pred)
    {
        Traverse_if(point, K, radius, AssignOp<T>(it), pred);
    }

    template <class ITERATOR, class PRED>
    void UnorderedExtract_if(const tPoint3<REAL>& point, REAL radius, ITERATOR it, const PRED& pred) const
    {
        KDTree* tree = const_cast<KDTree*>(this); // HACK
        tree->UnorderedTraverse_if(point, radius, AssignOp<T>(it), pred);
    }

    void render(float* cpoints, float* clines) const;
    void render(const tPoint3<REAL>& p, float* clines) const;

    void tree_depth(
        int& max_depth, 
        int& average_depth, 
        int& num_objects, 
        int& num_nodes, 
        int& num_leaves) const;

    // The implementation is here since there's a bug in VC6
    template <class IT>
    void Insert(IT first, IT last)
    {
        while (first < last)
        {
            //
            // Find the tree node that contains *first
            //
            tBox3<REAL> box;
            TreeNode* parent; // not used
            LeafNode* container = Find(_getpoint(*first), box, parent);

            //
            // find the rest of the objects that belong here
            //
            IT last_in_box = first+1;
            IT current = first+1;
            while (current < last)
            {
                if (box.contains(_getpoint(*current)))
                {
                    if (last_in_box != current)
                    {
                        // This does not work!!! to be implemented...
                        printf("KDTree::Insert Bug\n");
//                        std::iter_swap(last_in_box, current);
                    }
                    ++last_in_box;
                }
                ++current;
            }
            container->Insert(first, last_in_box);
            first = last_in_box;
        } 
    } // Insert

    //
    // Same as insert, but works for my piter iterator
    // adaptor.
    // A hack. Is there a clean solution?
    //
    template <class IT>
    void InsertSpecial(IT first, IT last)
    {
        typedef typename std::vector<IT> t_itlist;
        typedef typename t_itlist::iterator it2;

        t_itlist itlist(last-first);

        {
            it2 b = itlist.begin();
            for (; first != last; ++first, ++b)
            {
                *b = first;
            }
        }

        it2 ifirst = itlist.begin();
        it2 ilast = itlist.end();
        while (ifirst < ilast)
        {
            //
            // Find the tree node that contains *first
            //
            tBox3<REAL> box;
            TreeNode* parent; // not used
            LeafNode* container = Find(_getpoint(**ifirst), box, parent);

            //
            // find the rest of the objects that belong here
            //
            it2 last_in_box = ifirst+1;
            it2 current = ifirst+1;
            while (current < ilast)
            {
                if (box.contains(_getpoint(**current)))
                {
                    if (last_in_box != current)
                    {
                        std::iter_swap(last_in_box, current);
                    }
                    ++last_in_box;
                }
                ++current;
            }

            //
            // Documentation:
            // container->Insert(ifirst, last_in_box);
            // ifirst = last_in_box;
            //
            int P = container->Resize(last_in_box - ifirst);
            for (; ifirst != last_in_box; ++ifirst, ++P)
            {
                assert(**ifirst);
                (*container)[P] = **ifirst;
            }
        } 
    } // InsertSpecial

protected:
/*--------------------- Misc -----------------*/
    struct TruePred
    {
        bool operator()(const T&) const { return true; }
    };

public: // was protected: because InsertSpecial uses Node.
#if PROFILE_UNORDERED_QUERY==1
int profile_unorered_query_visited_nodes;
#endif

    /*--------------------- tree data structures -----------------*/

    struct Node
    {
        Node() {}
        virtual ~Node() {}
        typedef enum {leaf, tree, single} node_type; // Types of nodes

        virtual node_type Type() const=0;
        virtual void print_type() const=0;
    };

    struct TreeNode : public Node
    {
	typedef typename Node::node_type node_type;
        ~TreeNode() { delete l; delete r; }
        node_type Type() const { return Node::tree; }
        void print_type() const { printf("TreeNode"); }

        // Which axis is cut by this node
        // Numbering complies with VX, VY, VZ defined in algebra3.h
        enum t_axis {XAxis=0, YAxis, ZAxis}; 

        REAL cut_point;       // where along the axis is this node cut.
                                // Absolute coordinates.

        Node *l, *r;            // Left and right sons
        t_axis axis;
    }; // TreeNode

    struct LeafNode : public Node
    {
	typedef typename Node::node_type node_type;
        LeafNode() {}
        node_type Type() const { return Node::leaf; }
        void print_type() const { printf("LeafNode"); }
        
        void Insert(T v) { objects.push_back(v); }

        template<class IT>
        void Insert(IT first, IT last)
        {
            objects.insert(objects.end(), first, last);
        }

		// JS
		void Remove(T v) {
			if (objects.size() == 1) {
				if (objects[0] != v)
					printf("KDTree::Remove not found\n");
				else
					objects.resize(0);

			} else {
				for (unsigned i=0; i<objects.size(); i++) {
					if (objects[i] == v) {
						objects[i] = objects.back();
						objects.pop_back();
						return;
					}
				}

				printf("KDTree::Remove not found\n");
			}
		}

		bool empty() const { return objects.size() == 0; }

        //
        // Reserve space for additional K items, return previous # of elements
        //
        int Resize(int K)
        {
            int n = objects.size();
            objects.resize(n+K);
            return n;
        }

//        typedef tightvec<T> t_objects_list;
        typedef std::vector<T> t_objects_list;
        t_objects_list objects; // List of objects in the node

        typename t_objects_list::value_type& operator[](int n) {
                return objects[n];
        }

    }; // LeafNode

    struct SingleObject : public Node
    {
	typedef typename Node::node_type node_type;
        T* object;
        explicit SingleObject() : Node() {}
        SingleObject(T* r_object) : Node(), object(r_object) {}
        SingleObject(const SingleObject& rhs) : Node(), object(rhs.object) {}

        void set(T* r_object) { object = r_object; }
        node_type Type() const { return Node::single; }
        void print_type() const { printf("SingleObject"); }
    };
protected:
    /*--------------------- [ tree data structures ] -----------------*/

    /*--------------------------------------------------------------*/
    //
    // Search Algorithm:
    //
    // Traverse the points in the given neighborhood in increasing
    // distance from the give point.
    // Stop when reached maximal radius or maximal number of points.
    //
    // How it works:
    //   maintain a priority queue of cells / points
    //   at every iteration extract the closest object
    //   if it is a point, apply OP(p)
    //   otherwise expand the cell and insert it's contents into the queue.
    //   priority = distance.
    //
    //   Phase 1: Insert the relevant cells from the root to the cell of
    //            the center point
    //   Phase 2: Extract objects from the priority queue until finish.
    //


public:
    //
    // An object that holds information for traversing the
    // tree incrementally, i.e. return one item at a time
    // ordered by distance from the origin
    //
    class OrderedIncrementalTraverse
    {
    public:
        OrderedIncrementalTraverse(const KDTree& tree, const tPoint3<REAL>& origin) :
            _origin(origin),
            _tree(tree), 
            single_objects_buffer(100),
            allocated_single_objects_buffer(0)
        {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_ordered_init_pq);
#endif // __KDT_PROFILING==1
            InitPQ();
        }
      
        bool empty() const { return _pq.empty(); }

        const T& GetNext(REAL& squared_distance)
        {
#if __KDT_PROFILING==1
            HighResAutoStop _t(_timer_ordered_getnext);
#endif // __KDT_PROFILING==1
            while (!_pq.empty())
            {
                NeighborCell cell = _pq.top();
                _pq.pop();
#if DEBUG_TRAVERSE==1
                cell.print_type(); printf(": dist: %g\n", cell.GetDist());
                _pq.print();
#endif // DEBUG_TRAVERSE

                //
                // BUGBUG: I think it sould have been nicer if this switch 
                // would have been a call to virtual function, but one cannot
                // make a templated virtual function.
                // What is the solution?
                //
                switch (cell.GetNode()->Type())
                {
                case Node::single:
                {
                    SingleObject* single = static_cast<SingleObject*>(cell.GetNode());
                    squared_distance = cell.GetDist();
                    return *single->object;
                    break;
                }
                case Node::leaf:
                {
                    LeafNode* leaf = static_cast<LeafNode*>(cell.GetNode());
                    TraverseLeaf(leaf);
                    break;
                }
                case Node::tree:
                {
                    TreeNode* tree = static_cast<TreeNode*>(cell.GetNode());

                    REAL d = tree->cut_point - _origin[tree->axis];
                    Vector3 leftdistvector(cell.GetDistvector());
                    Vector3 rightdistvector(cell.GetDistvector());
                    if (d > 0)
                    {
                        // we're going left
                        rightdistvector[tree->axis] = -d;
                    }
                    else
                    {
                        leftdistvector[tree->axis] = d;
                    }

                    // BUGBUG: performance, can move the above inside the if
                    // and not compute one of the squared_length()!!
                    _pq.push( NeighborCell(tree->l, leftdistvector) );
                    _pq.push( NeighborCell(tree->r, rightdistvector) );
                    break;
                }
                } // Switch
            } // While
#ifndef NO_EXCEPTIONS
            throw _exception_no_more_values;
#endif
            assert(0);
//            throw CErr("No more values");
        }

    protected:
        //
        // A helper class: priority queue object.
        //
        struct NeighborCell
        {
        public:
            NeighborCell(Node* NODE,  const Vector3& dist) :
                node(NODE), distvector(dist)
            {
                assert(NODE);
//                assert(unsigned(NODE) < 0xf000000);
//                if (!owner) node.release();
                mindist = distvector.squared_length();
            }

            NeighborCell(Node* NODE, REAL dist) :
                node(NODE),  mindist(dist)
            {
//                if (!owner) node.release();
            }

            bool operator < (const NeighborCell& rhs) const
            {
#if FASTPQ2==0
                return mindist > rhs.mindist;
#else
                return mindist < rhs.mindist; // For std::set priority queue implementation.
#endif
            }

//            Node* GetNode() { return node.get(); }
            Node* GetNode() { return node;}
            const Vector3& GetDistvector() const { return distvector; }
            REAL GetDist() const { return mindist; }
            void print_type() const { node->print_type(); }

        private:
//            auto_ptr2<Node> node;
            Node* node;
            Vector3 distvector;
            REAL mindist; // minimal distance to cell on each of the axes

        };

//        typedef std::priority_queue<NeighborCell> t_pq;
#if FASTPQ2==1
        typedef fast_pq2<NeighborCell> t_pq;
#else
        typedef fast_pq<NeighborCell> t_pq;
#endif

        const tPoint3<REAL>& _origin;
        t_pq _pq;
        const KDTree& _tree;

        //
        // Initialize the priority queue as describe bellow under
        //  Search Algorithm:
        //
        //          Traverse the tree till we reach the cell
        //          the contains the point.
        //          At every node, if we go left, insert the right
        //             and if we got right, insert the left child.
        //
        //
        void InitPQ()
        {
            Node* cell = _tree.root;
            {
#if __KDT_PROFILING==1
                HighResAutoStop _t(_timer_ordered_init_pq2);
#endif // __KDT_PROFILING==1

                while (cell->Type() != Node::leaf)
                {
                    TreeNode* tree = static_cast<TreeNode*>(cell);

                    REAL d = tree->cut_point - _origin[tree->axis];

                    if (d > 0)
                    {
                        // we're going left
                        Vector3 dist(0.0);
                        dist[tree->axis] = -d;
                        _pq.push_only( NeighborCell(tree->r, dist ) );
                        cell = tree->l;
                    }
                    else
                    {
                        // going right
                        Vector3 dist(0.0);
                        dist[tree->axis] = d;
                        _pq.push_only( NeighborCell(tree->l, dist) );
                        cell = tree->r;
                    }
                } // While
            }

            Vector3 dist(0.0);
            {
				LeafNode* leaf = (LeafNode*)cell;
				
                if (!leaf->empty()) _pq.push_only( NeighborCell(cell, dist) );
            }
            {
                _pq.remake_heap();
            }
        }  // InitPQ

        void TraverseLeaf (LeafNode* leaf)
        {
            typename LeafNode::t_objects_list::iterator f = leaf->objects.begin();
            typename LeafNode::t_objects_list::iterator l = leaf->objects.end();
            for (; f != l; ++f)
            {
                REAL d = (_tree._getpoint(*f) - _origin).squared_length();

                SingleObject* single = allocate_single_object();
                single->set(&(*f));

                _pq.push(NeighborCell(single, d));
            }
        } // TraverseLeaf

        //
        // A test: replace allocating many SingleObject while traversing
        // This reduces traversal time by over 50 percent.
        // 
        // single_objects_buffer must be a container that does not reallocates!!!
        //
        std::deque<SingleObject> single_objects_buffer;
        unsigned allocated_single_objects_buffer;
        SingleObject* allocate_single_object()
        {
            if (allocated_single_objects_buffer >= single_objects_buffer.size()-1)
            {
                SingleObject dummy;
                single_objects_buffer.insert(single_objects_buffer.end(), 100, dummy);
            }

            SingleObject* r = &(single_objects_buffer[allocated_single_objects_buffer]);
            ++allocated_single_objects_buffer;
            return r;
        }
    }; // OrderedIncrementalTraverse

    OrderedIncrementalTraverse* NewOrderedIncrementalTraverse(const tPoint3<REAL>& origin) const
    {
        return new OrderedIncrementalTraverse(*this, origin);
    }

    /*-------------------- [ search algorithm ] ------------------------------*/
    //
    // Walk the tree - calls op() for each point
    // by the order of the tree
    //
    template <class FUNC> 
    void Walk(const FUNC& op)
    {
        if (root != 0) Walk(root, op);
    }
    template <class FUNC> void Walk(Node* parent, const FUNC& op);



    //
    // Search for the leaf node that contains p
    //
    // Return:
    //   bbox - for the leaf node containing p
    //   parent
    //
    LeafNode* Find(const tPoint3<REAL>& p, TreeNode*& parent) const;
    LeafNode* Find(const tPoint3<REAL>& p, tBox3<REAL>& bbox, TreeNode*& parent) const;

    void MakeTree(TreeNode* root);

protected: // JS moved

    //
    // Split a leaf node, returning a pointer to the
    // parent of the node that holds a new subtree
    //
    Node* Split (
        LeafNode* cell         // Input cell to split
        );

    /*-------------------- search algorithm ------------------------------*/
public:
    //
    // Don't know how to to declare member template function outside
    // of the class... Now I know, this is a vc bug...
    //
    template <class FUNC>
    void Traverse(
        const tPoint3<REAL>& point, 
        size_type K,
        REAL radius, // Maximal radius
        const FUNC& op) const
    {
        OrderedIncrementalTraverse* itrav = NewOrderedIncrementalTraverse(point);

        size_type n_visited = 0;   // Number of visited points so far.
        REAL maximal_radius = 0; // Maximal radius encountered so far

        REAL radius2 = radius * radius;

        while (!itrav->empty() && (n_visited < K) && (maximal_radius <= radius2) )
        {
            // BUGBUG: if radius is the limit may call with one more...
            op(itrav->GetNext(maximal_radius));
            ++n_visited;
        } // While

        delete itrav;
    } // Traverse

    //
    // Traverse the first K, or in given radius if PRED(point) holds
    // NOTE: if PRED(point) does not hold, points is not counted as
    //       one of the K points.
    //
    template <class FUNC, class PRED>
    void Traverse_if(
        const tPoint3<REAL>& center, 
        size_type K,
        REAL radius, // Maximal radius
        const FUNC& op,
        const PRED& pred)
    {
        aptr<OrderedIncrementalTraverse> itrav(NewOrderedIncrementalTraverse(center));

        size_type n_visited = 0;   // Number of visited points so far.
        REAL maximal_radius = 0; // Maximal radius encountered so far
        REAL radius2 = radius * radius;

        try
        {
            while (!itrav->empty() && (n_visited < K) && (maximal_radius <= radius2) )
            {
                const T& point = itrav->GetNext(maximal_radius);
                if (pred(point))
                {
                    op(point);
                    ++n_visited;
                }
            } // While
        }
        catch (kdtree_exception_no_more_values& no_more)
        {
            //(void)no_more;
        }
    } // Traverse_if

    template <class FUNC> void UnorderedTraverse(const tPoint3<REAL>& point, REAL radius, const FUNC& op)
    {
#if PROFILE_UNORDERED_QUERY==1
        profile_unorered_query_visited_nodes=0;
#endif
        TruePred alwaystrue;
        UnorderedTraverse_if(root, point, radius*radius, Vector3(0.0), op, alwaystrue);
    }

    template <class FUNC, class PRED> 
    void UnorderedTraverse_if(const tPoint3<REAL>& point, REAL radius, const FUNC& op, const PRED& pred)
    {
#if PROFILE_UNORDERED_QUERY==1
        profile_unorered_query_visited_nodes=0;
#endif
        UnorderedTraverse_if(root, point, radius*radius, Vector3(0.0), op, pred);
    }

private:

#if 0
    struct uoti_item
    {
        uoti_item() {}
        uoti_item(Node* node_) : node(node_), fix_item(false) {}
        uoti_item(REAL saved_coord_, TreeNode::t_axis axis_) : saved_coord(saved_coord_), axis(axis_), fix_item(true) {}
        union
        {
            REAL saved_coord;
            struct
            {
                TreeNode::t_axis axis;
                Node* node;
            };
        };
        bool fix_item; // true - means this is a return, and dist[axis] = saved_coord should be executed
                       // false - means this is a standard node to traverse
    };

    // Iterative version - works slow, need to check why...
    template <class FUNC, class PRED> 
    void UnorderedTraverse_if(Node* node, const tPoint3<REAL>& point, REAL radius2, Vector3& dist, const FUNC& op, const PRED& pred)
    {
        // The stack
        std::vector<uoti_item> l_stack;
        l_stack.reserve(100);
        l_stack.push_back(uoti_item(node));
        while (!l_stack.empty())
        {
#if PROFILE_UNORDERED_QUERY==1
    ++profile_unorered_query_visited_nodes;
#endif
            uoti_item& pnode = l_stack.back(); // HACK
            l_stack.pop_back();
            switch (pnode.fix_item)
            {
            case true:
                dist[pnode.axis] = pnode.saved_coord;
                break;
            case false:
                if (pnode.node->Type() == Node::leaf)
                {
                    LeafNode* leaf = static_cast<LeafNode*>(pnode.node);
                    typename LeafNode::t_objects_list::iterator f = leaf->objects.begin();
                    typename LeafNode::t_objects_list::iterator l = leaf->objects.end();
                    for (; f != l; ++f)
                    {
                        REAL d = (_getpoint(*f) - point).squared_length();

                        if ( (d <= radius2) && pred(*f) )
                        {
                            op(*f);
                        }
                    }
                }
                else
                {
                    TreeNode* tree = static_cast<TreeNode*>(pnode.node);

                    REAL d = tree->cut_point - point[tree->axis];

                    if (d > 0)
                    {
                        l_stack.push_back(uoti_item(tree->l));
                        l_stack.push_back(uoti_item(dist[tree->axis], tree->axis));
                        dist[tree->axis] = -d;
                        if (dist.squared_length() <= radius2)
                        {
                            l_stack.push_back(uoti_item(tree->r));
                        }
                    }
                    else
                    {
                        l_stack.push_back(uoti_item(tree->r));
                        l_stack.push_back(uoti_item(dist[tree->axis], tree->axis));
                        dist[tree->axis] = d;
                        if (dist.squared_length() <= radius2)
                        {
                            l_stack.push_back(uoti_item(tree->l));
                        }
                    }
                }
                break;
            } // swith pnode.fix_type
        } // while
    } // UnorderedTraverse_if (iterative version

#else // recursive version
    //
    // dist is the distance to the near this cell in each axis
    // so the real distance is dist.length().
    //
    template <class FUNC, class PRED> 
    void UnorderedTraverse_if(Node* node, const tPoint3<REAL>& point, REAL radius2, Vector3& dist, const FUNC& op, const PRED& pred)
    {
#if PROFILE_UNORDERED_QUERY==1
        ++profile_unorered_query_visited_nodes;
#endif
        if (node->Type() == Node::leaf)
        {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            typename LeafNode::t_objects_list::iterator f = leaf->objects.begin();
            typename LeafNode::t_objects_list::iterator l = leaf->objects.end();
            for (; f != l; ++f)
            {
                REAL d = (_getpoint(*f) - point).squared_length();

                if ( (d <= radius2) && pred(*f) )
                {
                    op(*f);
                }
            }
        }
        else
        {
            TreeNode* tree = static_cast<TreeNode*>(node);


            REAL d = tree->cut_point - point[tree->axis];

            if (d > 0)
            {
                // we're going left
                UnorderedTraverse_if(tree->l, point, radius2, dist, op, pred);
                REAL saved_coord = dist[tree->axis];
                dist[tree->axis] = -d;
                if (dist.squared_length() <= radius2)
                {
                    UnorderedTraverse_if(tree->r, point, radius2, dist, op, pred);
                }
                dist[tree->axis] = saved_coord;
            }
            else
            {
                UnorderedTraverse_if(tree->r, point, radius2, dist, op, pred);
                REAL saved_coord = dist[tree->axis];
                dist[tree->axis] = d;
                if (dist.squared_length() <= radius2)
                {
                    UnorderedTraverse_if(tree->l, point, radius2, dist, op, pred);
                }
                dist[tree->axis] = saved_coord;
            }

        }

#if 0
          // Same as above, but slower? for reference only 
            TreeNode* tree = static_cast<TreeNode*>(node);


            REAL d = tree->cut_point - point[tree->axis];

            Vector3 leftdistvector(dist);
            Vector3 rightdistvector(dist);

            if (d > 0)
            {
                // we're going left
                rightdistvector[tree->axis] = -d;
            }
            else
            {
                leftdistvector[tree->axis] = d;
            }

            if (leftdistvector.squared_length() <= radius2)
            {
                UnorderedTraverse_if(tree->l, point, radius2, leftdistvector, op, pred);
            }

            if (rightdistvector.squared_length() <= radius2)
            {
                UnorderedTraverse_if(tree->r, point, radius2, rightdistvector, op, pred);
            }
        }
#endif // 0
    } // UnorderedTraverse_if
#endif //0 recursive version 

    /*-------------------- [ search algorithm ] ------------------------------*/

private:

    void Print(Node* root, const tBox3<REAL>& bbox, int depth) const;
    void render(float* cpoints, float* clines, Node* arg_root, const tBox3<REAL>& arg_bbox) const;
    void render(const tPoint3<REAL>&p, float* clines, Node* arg_root, const tBox3<REAL>& arg_bbox) const;
    void tree_depth(Node* root_, int depth, int& max_depth, int& sum_depths, int& num_leaves, int& num_objects, int& num_nodes) const;

    size_type max_in_cell;            // Maximal # of objects in a cell
    tBox3<REAL> bbox;


public: // was protected: For OrderedIncrementalTraverse
    Node* root;
    GETPOINT _getpoint;

}; // KDTree

template<class T, class REAL, class GETPOINT>
inline KDTree<T, REAL, GETPOINT>::KDTree(size_type MAX_IN_CELL, const tBox3<REAL>& BBOX, const GETPOINT& getpoint) :
    max_in_cell(MAX_IN_CELL),
    bbox(BBOX),
    root(new LeafNode),
    _getpoint(getpoint)
{
}

template<class T, class REAL, class GETPOINT>
inline void KDTree<T, REAL, GETPOINT>::Insert(T object, bool split)
{
    TreeNode* parent;
    LeafNode* node = Find(_getpoint(object), parent);
    node->Insert(object);

    if (split)
    {
        if (parent)
        {
            if (parent->l == node) parent->l = Split(node);
            else parent->r = Split(node);
        }
        else
        {
            root = Split(node);
        }
    }
}

// JS
template<class T, class REAL, class GETPOINT>
inline void KDTree<T, REAL, GETPOINT>::Remove(T object)
{
    TreeNode* parent;
    LeafNode* node = Find(_getpoint(object), parent);
    node->Remove(object);
}


template<class T, class REAL, class GETPOINT>
inline T KDTree<T, REAL, GETPOINT>::FindMin(const tPoint3<REAL>& p) const
{
    SaveSingle sm;
    Traverse(p, 1, 1e30, sm);

    return sm.p;
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::MakeTree()
{
    if (root->Type() == Node::leaf)
    {
        root = Split(static_cast<LeafNode*>(root));
    }
    else
    {
        MakeTree(static_cast<TreeNode*>(root));
    }
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::MakeTree(TreeNode* arg_root)
{
    if (arg_root->l->Type() == Node::leaf)
    {
        arg_root->l = Split(static_cast<LeafNode*>(arg_root->l));
    }
    else
    {
        MakeTree(static_cast<TreeNode*>(arg_root->l));
    }

    if (arg_root->r->Type() == Node::leaf)
    {
        arg_root->r = Split(static_cast<LeafNode*>(arg_root->r));
    }
    else
    {
        MakeTree(static_cast<TreeNode*>(arg_root->r));
    }
}


template<class T, class REAL, class GETPOINT>
typename KDTree<T, REAL, GETPOINT>::LeafNode* KDTree<T, REAL, GETPOINT>::Find(const tPoint3<REAL>& p, typename KDTree<T, REAL, GETPOINT>::TreeNode*& parent) const
{
    if (root->Type() == Node::leaf)
    {
        parent = 0;
        return static_cast<LeafNode*>(root);
    }
    else
    {
        parent = static_cast<TreeNode*>(root);
        while (1)
        {
            if (parent->cut_point < p[parent->axis])
            {
                // we're going right
                if (parent->r->Type() == Node::leaf)
                {
                    return static_cast<LeafNode*>(parent->r);
                }
                parent = static_cast<TreeNode*>(parent->r);
            }
            else
            {
                // we're going left
                if (parent->l->Type() == Node::leaf)
                {
                    return static_cast<LeafNode*>(parent->l);
                }
                parent = static_cast<TreeNode*>(parent->l);
            }
        }
    }
}

template<class T, class REAL, class GETPOINT>
typename KDTree<T, REAL, GETPOINT>::LeafNode* KDTree<T, REAL, GETPOINT>::Find(const tPoint3<REAL>& p, tBox3<REAL>& leafbbox, typename KDTree<T, REAL, GETPOINT>::TreeNode*& parent) const
{
    leafbbox = bbox;
    if (root->Type() == Node::leaf)
    {
        parent = 0;
        return static_cast<LeafNode*>(root);
    }
    else
    {
        parent = static_cast<TreeNode*>(root);
        tBox3<REAL> lbbox, rbbox;
        while (1)
        {
            SplitBBox(leafbbox, parent->axis, parent->cut_point, lbbox, rbbox);
            if (parent->cut_point < p[parent->axis])
            {
                // we're going right
                leafbbox = rbbox;
                if (parent->r->Type() == Node::leaf)
                {
                    return static_cast<LeafNode*>(parent->r);
                }
                parent = static_cast<TreeNode*>(parent->r);
            }
            else
            {
                // we're going left
                leafbbox = lbbox;
                if (parent->l->Type() == Node::leaf)
                {
                    return static_cast<LeafNode*>(parent->l);
                }
                parent = static_cast<TreeNode*>(parent->l);
            }
        }
    }
}

template<class T, class REAL, class GETPOINT>
typename KDTree<T, REAL, GETPOINT>::Node* KDTree<T, REAL, GETPOINT>::Split (
    typename KDTree<T, REAL, GETPOINT>::LeafNode* cell  // Input cell to split
    )
{
    unsigned n_objects = cell->objects.size();
    if (n_objects <= max_in_cell) return cell;

        //
        // Compute a split point and axis with the following heuristic:
        // Choose three points at random
        // compute their bounding box
        // Take the split axis to be the maximum of the width,height,depth
        //   of the bounding box
        // Take the position of the split at the middle of the bounding box
        //
        typename TreeNode::t_axis axis;
        REAL cut_point;
        {
                unsigned idx1 = rand32() % n_objects;
                unsigned idx2 = rand32() % n_objects;
                unsigned idx3 = rand32() % n_objects;
                const tPoint3<REAL>& p1 = _getpoint(cell->objects[idx1]);
                const tPoint3<REAL>& p2 = _getpoint(cell->objects[idx2]);
                const tPoint3<REAL>& p3 = _getpoint(cell->objects[idx3]);

                tPoint3<REAL> min_point(
                        min3(p1[0],p2[0],p3[0]),
                        min3(p1[1],p2[1],p3[1]),
                        min3(p1[2],p2[2],p3[2]));
                tPoint3<REAL> max_point(
                        max3(p1[0],p2[0],p3[0]),
                        max3(p1[1],p2[1],p3[1]),
                        max3(p1[2],p2[2],p3[2]));

                tVector3<REAL> box_dimentions = max_point - min_point;

                int k;
                max3(box_dimentions[0], box_dimentions[1], box_dimentions[2], k);
//                REAL an_epsilon = box_dimentions[k]/100;
                axis = static_cast<typename TreeNode::t_axis>(k);
                cut_point = min_point[axis] + box_dimentions[axis]*0.51;
        }

    //
    // Make the split
    //
//    lprintf("Splitting at %d %g\n", axis, cut_point);

    LeafNode* rightnode = new LeafNode;

    typename LeafNode::t_objects_list::iterator last_left = cell->objects.end();
    typename LeafNode::t_objects_list::iterator f = cell->objects.begin();
    while (f < last_left)
    {
        if (_getpoint(*f)[axis] <= cut_point)
        {
            ++f;
        }
        else
        {
            --last_left;
            rightnode->objects.push_back(*f);
            *f = *last_left;
        }
    }

    if (last_left == cell->objects.end())
    {
        // we have many objects that are equal
        delete rightnode;
        return cell;
    }

    cell->objects.erase(last_left, cell->objects.end());
    compact(cell->objects);

    TreeNode* subroot = new TreeNode;
    subroot->cut_point = cut_point;
    subroot->axis = axis;
    subroot->l = Split(cell);
    subroot->r = Split(rightnode);
    return subroot;
}

    //
    // Walk the tree - calls op() for each point
    // by the order of the tree
    //
template<class T, class REAL, class GETPOINT> template<class FUNC>
void KDTree<T, REAL, GETPOINT>::Walk(Node* parent, const FUNC& op)
{
    if (parent->Type() == Node::tree)
    {
        TreeNode* node = static_cast<TreeNode*>(parent);
        Walk(node->l, op);
        Walk(node->r, op);
    }
    else // if type == leaf
    {
        LeafNode* leaf = static_cast<LeafNode*>(parent);
        int K = leaf->objects.size();
        for (int i = 0; i < K; ++i)
        {
            op((*leaf)[i]);
        }
    }
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::render(float* cpoints, float* clines) const
{
    glPushAttrib(GL_POLYGON_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    render(cpoints, clines, root, bbox);
    glPopAttrib();
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::render(float* cpoints, float* clines, Node* arg_root, const tBox3<REAL>& arg_bbox) const
{
    if (arg_root->Type() == Node::leaf)
    {
        {
            glColor4fv(clines);
            arg_bbox.render();
        }
    }
    else
    {
        TreeNode* tree = static_cast<TreeNode*>(arg_root);
        tBox3<REAL> lbbox, rbbox;
        SplitBBox(arg_bbox, tree->axis, tree->cut_point, lbbox, rbbox);
        render(cpoints, clines, tree->l, lbbox);
        render(cpoints, clines, tree->r, rbbox);
    }
}

/*
 * Render bounding boxes from root to the point, for debug
 */
template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::render(const tPoint3<REAL>& p, float* clines) const
{
    glPushAttrib(GL_POLYGON_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor4fv(clines);
    render(p, clines, root, bbox);
    glPopAttrib();
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::render(const tPoint3<REAL>& p, float* clines, Node* arg_root, const tBox3<REAL>& arg_bbox) const
{
    arg_bbox.render();
    if (arg_root->Type() != Node::leaf)
    {
        TreeNode* tree = static_cast<TreeNode*>(arg_root);
        tBox3<REAL> lbbox, rbbox;
        SplitBBox(arg_bbox, tree->axis, tree->cut_point, lbbox, rbbox);
        if (lbbox.contains(p))
        {
            render(p, clines, tree->l, lbbox);
        }
        else
        {
            render(p, clines, tree->r, rbbox);
        }
    }
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::tree_depth(int& max_depth, int& average_depth, int& num_objects, int& num_nodes, int& num_leaves) const
{
    int sum_depths=0;
    num_leaves=0;
    max_depth = 0;
    num_objects = 0;
    num_nodes = 0;
    tree_depth(root, 0, max_depth, sum_depths, num_leaves, num_objects, num_nodes);
//    printf("Tree has %d leaves\n", num_leaves);
    average_depth = sum_depths / num_leaves;
}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::tree_depth(Node* root_, int depth, int& max_depth, int& sum_depths, int& num_leaves, int& num_objects, int& num_nodes) const
{
    ++num_nodes;

    if (root_->Type() == Node::leaf)
    {
        LeafNode* leaf = static_cast<LeafNode*>(root_);
        max_depth = max2(max_depth, depth);
        sum_depths += depth;
        ++num_leaves;
        num_objects += leaf->objects.size();
    }
    else
    {
        TreeNode* tree = static_cast<TreeNode*>(root_);
        int dl = 0;
        int dr = 0;
        tree_depth(tree->l, depth+1, dl, sum_depths, num_leaves, num_objects, num_nodes);
        tree_depth(tree->r, depth+1, dr, sum_depths, num_leaves, num_objects, num_nodes);
        max_depth = max2(dl, dr);
    }
}

template<class T, class REAL, class GETPOINT>
inline void KDTree<T, REAL, GETPOINT>::PrintStatistics() const
{
  // JS - FIXME
#ifdef WIN32
    printf("kdt: ordered query init pq: "); _timer_ordered_init_pq.print() printf("\n");
    printf("kdt: ordered query get next: "); _timer_ordered_getnext.print() printf("\n");
#endif
}

//
// DEBUG!!!
//
inline void xPrint(unsigned x) { printf("%d", x);}
template <class REAL> inline void xPrint(const tPoint3<REAL>& p) { printf("[%g %g %g]\n", p[0], p[1], p[2]);}

template<class T, class REAL, class GETPOINT>
void KDTree<T, REAL, GETPOINT>::Print(Node* arg_root, const tBox3<REAL>& arg_bbox, int depth) const
{
    if (arg_root->Type() == Node::leaf)
    {
        LeafNode* leaf = static_cast<LeafNode*>(arg_root);
        
        typename LeafNode::t_objects_list::iterator i = leaf->objects.begin();
        for (; i != leaf->objects.end(); ++i)
        {
            printf("%*c", depth*2, ' ');
            xPrint(*i);
        }
    }
    else
    {
        TreeNode* tree = static_cast<TreeNode*>(arg_root);
        tBox3<REAL> lbbox, rbbox;
        SplitBBox(arg_bbox, tree->axis, tree->cut_point, lbbox, rbbox);
        printf("%*c%d(%f)\n", depth*2,'L', tree->axis, tree->cut_point);
        Print(tree->l, lbbox, depth+1);
        printf("%*c%d(%f)\n", depth*2,'R', tree->axis, tree->cut_point);
        Print(tree->r, rbbox, depth+1);
    }
}

/*---------------------- Helper functions --------------*/
template <class REAL>
inline void SplitBBox(
    const tBox3<REAL>& arg_bbox,       // Input bounding box
    int axis,
    REAL position,            // axis, position define where to split
    tBox3<REAL>& lbbox, 
    tBox3<REAL>& rbbox                 // Result is store in [lr]bbox
    )
{
    assert(position <= arg_bbox.max_point()[axis]);
    assert(position >= arg_bbox.min_point()[axis]);

    lbbox = rbbox = arg_bbox;

    lbbox.set_max(axis, position);
    rbbox.set_min(axis, position);
}
/*---------------------- [ Helper functions ] --------------*/





/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/

// kd tree that allows overlapping bounding boxes - each object only exists in the tree in one place


// DISTBOXCLASS	must provide two typenames and two functions:
// typename DISTBOXCLASS::Box3
// typename DISTBOXCLASS::Point3
// DISTBOXCLASS::bounding_box(const OBJECT &) const;
// DISTBOXCLASS::distance(const OBJECT &, const typename DISTBOXCLASS::Point3 &) const;


template <typename OBJECT, typename DISTBOXCLASS>
class BoxKDTree {

public:

#define BOXKDTREE_MAXLEAFSIZE	10

	BoxKDTree() {
		children[0] = children[1] = NULL;
	}
	~BoxKDTree() {
		if (children[0]) delete children[0];
		if (children[1]) delete children[1];
	}


	BoxKDTree(const std::vector<OBJECT> &iobjects, const DISTBOXCLASS &boxclass, int axis=0) {

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
			bbox = DISTBOXCLASS::Box3::make_union(bbox, boxclass.bounding_box(objects[i]));
		}

		Split(boxclass, axis);

	}


	// get all the objects who's bounding boxes intersect the given box
	void GetIntersectedBoxes(const DISTBOXCLASS &boxclass, const typename DISTBOXCLASS::Box3 &ibox, std::vector<OBJECT> &intersected) const {

		// check our bounding box
		if (ibox.classify_position(bbox) == DISTBOXCLASS::Box3::OUTSIDE) {
			return;
		}


		// check any leaf objects
		for (unsigned i=0; i<objects.size(); i++) {
			typename DISTBOXCLASS::Box3 obox = boxclass.bounding_box(objects[i]);

			if (ibox.classify_position(obox) != DISTBOXCLASS::Box3::OUTSIDE) {
				intersected.push_back(objects[i]);
			}
		}

		// try going into the children
		if (children[0])	children[0]->GetIntersectedBoxes(boxclass, ibox, intersected);
		if (children[1])	children[1]->GetIntersectedBoxes(boxclass, ibox, intersected);

	}

	// get all the objects who's bounding boxes intersect the given box
	void GetIntersectedBoxes(const DISTBOXCLASS &boxclass, const typename DISTBOXCLASS::Point3 &p, std::vector<OBJECT> &intersected) const {

		// check our bounding box
		if (!bbox.contains(p)) {
			return;
		}


		// check any leaf objects
		for (unsigned i=0; i<objects.size(); i++) {
			typename DISTBOXCLASS::Box3 obox = boxclass.bounding_box(objects[i]);

			if (bbox.contains(p)) {
				intersected.push_back(objects[i]);
			}
		}

		// try going into the children
		if (children[0])	children[0]->GetIntersectedBoxes(boxclass, p, intersected);
		if (children[1])	children[1]->GetIntersectedBoxes(boxclass, p, intersected);

	}


	void Insert(const DISTBOXCLASS &boxclass, const OBJECT &o, int axis=0) {

		typename DISTBOXCLASS::Box3 obox = boxclass.bounding_box(o);

		// figure out which side we want to put it in
		int addside = -1;

		if (children[0] && children[1]) {
			// see which insertion would result in a smaller bounding box overlap

			typename DISTBOXCLASS::Box3 c0e = DISTBOXCLASS::Box3::make_union(children[0]->bbox, obox);
			typename DISTBOXCLASS::Box3 c1e = DISTBOXCLASS::Box3::make_union(children[1]->bbox, obox);

			bool intersect0 = c0e.classify_position(children[1]->bbox) != DISTBOXCLASS::Box3::OUTSIDE;
			bool intersect1 = c1e.classify_position(children[0]->bbox) != DISTBOXCLASS::Box3::OUTSIDE;

			if (intersect0 && !intersect1) {
				addside = 1;
			} else if (!intersect0 && intersect1) {
				addside = 0;
			} else if (intersect0 && intersect1) {
				// figure out which way causes the smallest overlap

			    typename DISTBOXCLASS::Point3
				t1(std::max(c0e.x_min(),children[1]->bbox.x_min()),
				   std::max(c0e.y_min(),children[1]->bbox.y_min()),
				   std::max(c0e.z_min(),children[1]->bbox.z_min())),
				t2(std::min(c0e.x_max(),children[1]->bbox.x_max()),
				   std::min(c0e.y_max(),children[1]->bbox.y_max()),
				   std::min(c0e.z_max(),children[1]->bbox.z_max()));
			    typename DISTBOXCLASS::Box3 ibox0(t1, t2);
			    typename DISTBOXCLASS::Point3
				t3(std::max(c1e.x_min(),children[0]->bbox.x_min()),
				   std::max(c1e.y_min(),children[0]->bbox.y_min()),
				   std::max(c1e.z_min(),children[0]->bbox.z_min())),
				t4(std::min(c1e.x_max(),children[0]->bbox.x_max()),
				   std::min(c1e.y_max(),children[0]->bbox.y_max()),
				   std::min(c1e.z_max(),children[0]->bbox.z_max()));
			    typename DISTBOXCLASS::Box3 ibox1(t3, t4);

				if (ibox0.x_length()*ibox0.y_length()*ibox0.z_length() < ibox1.x_length()*ibox1.y_length()*ibox1.z_length())
					addside = 0;
				else
					addside = 1;

			} else {
				// adding to neither would cause an intersection - add to the one that increases volume the least
				if (c0e.x_length()*c0e.y_length()*c0e.z_length() < c1e.x_length()*c1e.y_length()*c1e.z_length())
					addside = 0;
				else
					addside = 1;
			}



		} else if (children[0] && !children[1]) {
			addside = 0;
		} else if (!children[0] && children[1]) {
			addside = 1;
		}


		
		// expand our own bounding box
		bbox = (addside==-1 && objects.size()==0) ? obox : DISTBOXCLASS::Box3::make_union(bbox, obox);

		if (addside == -1) {
			objects.push_back(o);
			Split(boxclass, axis);
		} else {
			children[addside]->Insert(boxclass, o, (axis+1)%3);
		}

	}



	bool Remove(const DISTBOXCLASS &boxclass, const OBJECT &o) {

		if (bbox.classify_position(boxclass.bounding_box(o)) == DISTBOXCLASS::Box3::OUTSIDE)
			return false;


		// first check in the list of objects at this node
		for (unsigned i=0; i<objects.size(); i++) {
			if (o == objects[i]) {

				// remove the object from the list
				if (i != objects.size()-1) {
					objects[i] = objects.back();
				}
				objects.pop_back();

				// recompute the bounding box
				if (objects.size() > 0) {
			        bbox = boxclass.bounding_box(objects[0]);
					for (unsigned i=1; i<objects.size(); i++) {
						bbox = DISTBOXCLASS::Box3::make_union(bbox, boxclass.bounding_box(objects[i]));
					}
				} else {
				    typename DISTBOXCLASS::Point3
					u(-1,-1,-1),
					d(1,1,1);
				    typename DISTBOXCLASS::Box3
					b(u, d);
				    bbox = b;//DISTBOXCLASS::Box3(DISTBOXCLASS::Point3(-1,-1,-1), DISTBOXCLASS::Point3(1,1,1));
				}

				return true;
			}
		}


		// if we got here, we didn't find a match is the object list - check the children

		for (int c=0; c<2; c++) {
			if (children[c] && children[c]->Remove(boxclass, o)) {
				int dangle = children[c]->dangling();
				if (dangle != -1) {

					// the child we removed from now has no leaf objects and a single child - prune it from the tree
					BoxKDTree *child = children[c];
					children[c] = child->children[dangle];
					child->children[dangle] = NULL;
					delete child;

				} else if (children[c]->empty()) {
					// the child is now completely empty
					delete children[c];
					children[c] = NULL;
				}

				return true;
			}
		}

		// didn't find it anywhere!
		return false;
	}




	class OrderedTraverse {

	public:
		OrderedTraverse(BoxKDTree &root, const typename DISTBOXCLASS::Point3 &p, const DISTBOXCLASS &_distclass) : distclass(_distclass) {
			point = p;
			pq.push(std::pair<double, std::pair<bool,void*> >(-root.bbox.distance(p), std::pair<bool,void*>(false,&root)));
		}

		typename DISTBOXCLASS::Box3::value_type next(OBJECT &o) {
			
			while (!pq.empty()) {

				std::pair<double, std::pair<bool,void*> > top = pq.top();
				pq.pop();

				if (top.second.first) {
					// it's a pointer to a triangle
					o = *(OBJECT*)top.second.second;
					return -top.first; // we negate the disance so the smallest is at the top of the heap
				} else {
					// otherwise it's a kdtree node index
					BoxKDTree *node = (BoxKDTree*)top.second.second;

					// add any objects that we might have in this node
					for (unsigned i=0; i<node->objects.size(); i++) {
						pq.push(std::pair<double, std::pair<bool,void*> >(-distclass.distance(node->objects[i], point), std::pair<bool,void*>(true,&node->objects[i])));
					}

					// add the children if they exist
					if (node->children[0]) {
						pq.push(std::pair<double, std::pair<bool,void*> >(-node->children[0]->bbox.distance(point), std::pair<bool,void*>(false,node->children[0])));
					}

					if (node->children[1]) {
						pq.push(std::pair<double, std::pair<bool,void*> >(-node->children[1]->bbox.distance(point), std::pair<bool,void*>(false,node->children[1])));
					}

				}


			}

			return -1;

		}

	private:
		fast_pq< std::pair<double, std::pair<bool,void*> > >  pq;
		typename DISTBOXCLASS::Point3 point;

		const DISTBOXCLASS &distclass;
	};


private:
	void Split(const DISTBOXCLASS &boxclass, int axis=-1) {


		// check if we should stop splitting
		if (objects.size() <= BOXKDTREE_MAXLEAFSIZE) {
			compact(objects);
		} else {

			// if we're not told what axis to use, use the biggest
			if (axis == -1) {
				typename DISTBOXCLASS::Box3::value_type xl, yl, zl;
				xl = bbox.x_length();
				yl = bbox.y_length();
				zl = bbox.z_length();

				if (xl>yl && xl>zl) axis = 0;
				else if (yl>xl && yl>zl) axis = 1;
				else axis = 2;
			}


			// split the list by the axis
			std::vector<OBJECT> cobjects[2];

			if (objects.size() < 500) {
				std::vector< std::pair<double,OBJECT> > sorter(objects.size());

				for (unsigned i=0; i<objects.size(); i++) {
					typename DISTBOXCLASS::Box3 obox = boxclass.bounding_box(objects[i]);
					sorter[i] = std::pair<double,OBJECT>((double)obox.centroid()[axis], objects[i]);
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
					typename DISTBOXCLASS::Box3 obox = boxclass.bounding_box(objects[i]);

					if (obox.centroid()[axis] < bbox.centroid()[axis]) {
						cobjects[0].push_back(objects[i]);
					} else {
						cobjects[1].push_back(objects[i]);
					}
				}
			}


			if ((cobjects[0].size() != 0) && (cobjects[1].size() != 0)) {

				// actually have to split
				objects.clear();
				compact(objects);

				assert(!children[0] && !children[1]);

				children[0] = new BoxKDTree(cobjects[0], boxclass, (axis+1)%3);
				children[1] = new BoxKDTree(cobjects[1], boxclass, (axis+1)%3);

			}
		}
	}


	bool empty() {
		return (!children[0] && !children[1] && !objects.size());
	}

	int dangling() {
		if (!objects.size()) {
			if (children[0] && !children[1])	return 0;
			if (!children[0] && children[1])	return 1;
		}
		return -1;
	}

	// leaf node
	std::vector<OBJECT> objects;

	// internal node
	BoxKDTree* children[2];

	typename DISTBOXCLASS::Box3 bbox;

};









GTB_END_NAMESPACE

#endif // __KDTREE_H
