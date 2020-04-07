
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




// the front is just a list of vertices.  by arbitrary convention, the area still to be triangulated
// is to the left of the front (looking at the next vertex in the list, with surface normal up).



#ifndef _FRONT_H
#define _FRONT_H

#include <list>
#include "common.h"

class Front;



// these three are mutually exclusive
#define PRIORITY_GROW_SEED_EDGE     0x01000000
#define PRIORITY_GROW_EDGE		    0x02000000
#define PRIORITY_CONNECT		    0x04000000

// can be combined with any of these
#define PRIORITY_OWA			    0x08000000		// outside working area
#define PRIORITY_OWB			    0x10000000		// outside working block

#define PRIORITY_FAILSAFE		    0x20000000
#define PRIORITY_FAILSAFE_WILLFAIL	0x40000000
#define PRIORITY_FAILSAFE_ANY		(PRIORITY_FAILSAFE | PRIORITY_FAILSAFE_WILLFAIL)

// trumps all
#define PRIORITY_CONNECTOR          0x80000000


#define FRONT_FLAG_BOUNDARY		1   // edge's projection when onto boundary (hack for when growth goes outside the domain)
#define FRONT_FLAG_FEATURE		2   // sharp edges or boundaries
#define FRONT_FLAG_CONNECTOR    4   // pseudo-edge connecting boundary verts that will be filled in by boundaries in future blocks
#define FRONT_FLAG_EDGEGROW     8   // was the triangle behind this edge the result of an edge grow (vs connection)
#define FRONT_FLAG_SEED         16  // an edge that should be grown asap (seed points, boundaries, features)


class ProjectionResult {
    public:
    Point3 position;
    Vector3 normal;
	real_type step;
    volatile int result;   // success/boundary/fail/notfinished...
};

class FrontElement {
    public:

    typedef std::pair<unsigned int,real_type> priority_type;


    FrontElement() { }
    FrontElement(const Point3 &p, const Vector3 &n, int vi, real_type step)
	: position(p),
	  normal(n),
	  vertindex(vi),
	  priority(0x7fffffff,1e34),
	  max_step(step),
	  flags(0)
    {
    }


    // info about the vertex
    Point3 position;
    Vector3 normal;

    // the index for the output file (triangles share vertices)
    int vertindex;
    int vertindex_l;

    // the priority of the edge starting at this vertex
    priority_type priority;
	
    // how far we can step from this spot
    real_type max_step;

    int flags;

    ProjectionResult proj_res;


    // keep track of the position in the heap so we can update it's priority
    int heap_position;

    // we need to know which front we're actually a part of
    Front *front;
};



class Front {

    public:

    typedef std::list<FrontElement> fel;	// front element list
    typedef fel::iterator feli;		// front element list iterator, basically just pointers to FrontElement's,
                                        // nice because they don't get invalidated when shuffeling them around in the lists


    Front();

    // add an element to the end of the list - useful for creating the inital fronts
    feli AddElement(FrontElement &e);



    // split the front at verts i1 and i2, return a pointer to the new front created
    static Front* Split(feli i1, feli i2, feli &n1, feli &n2);

    // merge two fronts, i1's front grows, i2's front will need to be deleted
    static void Merge(feli i1, feli i2, feli &n1, feli &n2);


    // get the next/prev element on the front
    feli FirstElement();
    static feli NextElement(feli i);
    static feli PrevElement(feli i);
    static void RemoveElement(feli i);
    static feli InsertElement(feli next, FrontElement &e);

    bool empty() { return elements.empty(); }
    void verify();

    private:
    // don't allow copies
    Front(const Front &rhs) {
	cerr<<"trying to copy a front"<<endl;
	BREAK;
    }
    Front& operator=(const Front &rhs) {
	cerr<<"trying to copy a front"<<endl;
	BREAK;
	return *this;
    }


    fel elements;

};




#endif





