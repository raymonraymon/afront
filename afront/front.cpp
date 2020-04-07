
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



#include "common.h"
#include "front.h"


Front::Front() {
}

Front::feli Front::AddElement(FrontElement &e) {
    e.front = this;
    elements.push_back(e);
    return (--elements.end());
}



// split the front at verts i1 and i2, return a pointer to the new front created
Front* Front::Split(feli i1, feli i2, feli &n1, feli &n2) 
{
    if (i1->front != i2->front) {
	cerr<<"splitting different fronts!"<<endl;
	BREAK;
    }

    Front *f1 = i1->front;
    Front *f2 = new Front();

    feli cutstart, cutend;

    cutstart = cutend = i2;

    while (cutend != f1->elements.end() && cutend !=i1) {
	cutend->front = f2;
	++cutend; 
    }


    if (cutend == f1->elements.end()) {
	f2->elements.splice(f2->elements.end(), f1->elements, cutstart, cutend);

	// got to the end before finding i1, have to start again at the beginning
	cutstart = cutend = f1->elements.begin();

	while (cutend != f1->elements.end() && cutend !=i1) {
	    cutend->front = f2;
	    ++cutend;
	}

	if (cutend == f1->elements.end()) {
	    cerr<<"cutend at end!"<<endl;
	    BREAK;
	}

	f2->elements.splice(f2->elements.end(), f1->elements, cutstart, cutend);
    } else {
	f2->elements.splice(f2->elements.end(), f1->elements, cutstart, cutend);
    }


    // i1 is still in f1, i2 is now in f2

    // add a copy of i1 to f2 just before i2
    n1 = i2;
    f2->elements.insert(n1, *i1);
    --n1;
    n1->front = f2;

    // add a copy of i2 to f1 just before i1
    n2 = i1;
    f1->elements.insert(n2, *i2);
    --n2;
    n2->front = f1;


    return f2;
}



// merge two fronts, i1's front grows, i2's front will need to be deleted
void Front::Merge(feli i1, feli i2, feli &n1, feli &n2) {

    if (i1->front == i2->front) {
	cerr<<"Merging a single front"<<endl;
	BREAK;
    }

    Front *f1 = i1->front;
    Front *f2 = i2->front;

    // reset the fronts that f2 is pointing to
    for (feli i=f2->elements.begin(); i!=f2->elements.end(); ++i) {
	i->front = f1;
    }

    // copy from i2 to the end of f2 to just before i1
    f1->elements.splice(i1, f2->elements, i2, f2->elements.end());

    // copy the rest of the list
    f1->elements.splice(i1, f2->elements, f2->elements.begin(), f2->elements.end());


    // add a copy of i1 to f2 just before i2
    n1 = i2;
    f2->elements.insert(n1, *i1);
    --n1;

    // add a copy of i2 to f1 just before i1
    n2 = i1;
    f1->elements.insert(n2, *i2);
    --n2;
}






// get the next/prev element on the front
Front::feli Front::NextElement(feli i) {

    Front *f = i->front;

    ++i;

    // if we went off the end, go back to the beginning
    if (i == f->elements.end())
	return f->elements.begin();

    return i;
}

Front::feli Front::PrevElement(feli i) {

    Front *f = i->front;

    // if we're at the beginning, go back to the end
    if (i == f->elements.begin())
	return (--f->elements.end());

    return --i;
}

Front::feli Front::FirstElement() {
    return elements.begin();
}


Front::feli Front::InsertElement(feli next, FrontElement &e) {
    Front *f = next->front;
    e.front = f;
    f->elements.insert(next, e);
    return --next;
}

void Front::RemoveElement(feli i) {
    Front *f = i->front;
    f->elements.erase(i);
}




void Front::verify() {

    if (empty())
	return;

    feli e = FirstElement();

    do {

	feli next = NextElement(e);
	feli next2 = NextElement(next);

	if (e->vertindex == next->vertindex) {
	    cerr<<"next vertindex the same!"<<endl;
	    BREAK;
	}

	if (e->flags & FRONT_FLAG_CONNECTOR && next->flags & FRONT_FLAG_CONNECTOR) {
		cerr<<"two connectors in a row!"<<endl;
		BREAK;
	}

	/*
	if (e->vertindex == next2->vertindex) {
	    cerr<<"next2 vertindex the same!"<<endl;
	    BREAK;
	}
	*/

	real_type dist = Point3::distance(e->position, next->position);
	if (dist < e->max_step*0.1) {
	    cerr<<"close front vertex locations!"<<dist<<endl;

	    dbgClear();
	    
	    DbgPoints::add(e->position, 1,0,0);
	    DbgPoints::add(next->position, 1,0,0);

	    redrawAndWait(' ', true);
	}
	    

	e = next;

    } while (e != FirstElement());

}

