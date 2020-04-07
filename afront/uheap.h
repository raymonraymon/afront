
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



#ifndef _UHEAP_H
#define _UHEAP_H


#include <vector>
#include "common.h"
#include <cassert>

// updatable heap

// heap class that calls a function heap_position(T*, int) when an element gets moved to a
// new position in the heap.  this allows you to track where an element is, and update the priority
// for just that element, which can be done in log time instead of remaking the whole heap in linear time


// maintains the heap property of each node's children being less than the node


// requires the following two global functions defined:
// bool operator<(const T &, const T &);
// void heap_position(T &, int);



template <typename T>
class UHeap {

public:

	UHeap() { }
	~UHeap() { }

	// put in at the end, and move it up the heap
	void push(T &x) {
		int level = heap_level(x);
		int position = heap[level].size();

		if (level==1)
			position = -position-1;

		heap_position(x, position);
		heap[level].push_back(x);
		update_position(position);
	}



	// move the back to the front, and push it back down the heap
	void pop() {
		if (!heap[0].empty())
			remove(0);
		else {
			remove(-1);
		}
	}


	const T& top() {
		if (!heap[0].empty())
			return heap[0].front();
		else {
			return heap[1].front();
		}
	}



	void remove(const int i) {

		int level = 0;
		int li = i;
		if (li < 0) {
			li = -li-1;
			level = 1;
		}

		std::vector<T> &lheap = heap[level];

		heap_position(lheap[li], 0xFFFFFFFF);

		// swap the last element into this location
		if (li != (int)lheap.size()-1) {
			lheap[li] = lheap.back();
			heap_position(lheap[li], i);
		}
		lheap.pop_back();

		// fix the heap
		if ((lheap.size() > 0) && (li != lheap.size()))
			update_position(i);
	}


	void update_position(const int i) {
		int cur_level = 0;
		int li = i;
		if (li < 0) {
			li = -li-1;
			cur_level = 1;
		}

		int want_level = heap_level(heap[cur_level][li]);
		if (cur_level == want_level) {
			// update in place
			update_position(cur_level, li);
		} else {
			// moving from one level to the other - just remove it and re-add it
			T x = heap[cur_level][li];
			remove(i);
			push(x);
		}
	}


	bool empty() const {
		return heap[0].empty() && heap[1].empty();
	}

    int size()
    {
        return heap[0].size() + heap[1].size();
    }

	// void verify(int i=0) const {

	// 	if (i<(int)heap.size()) {

    //         assert(heap[i]->heap_position == i);

	// 		int p=parent(i);
	// 		if (p>=0) {
	// 		    if (heap[p] < heap[i]) {
	// 			cerr<<"heap out of order"<<endl;
	// 			BREAK;
	// 		    }
	// 		}

	// 		int lc=left(i);
	// 		if (lc<(int)heap.size()) {
	// 			if (heap[i] < heap[lc]) {
	// 			cerr<<"heap out of order"<<endl;
	// 			BREAK;
	// 		    }
	// 			verify(lc);
	// 		}

	// 		int rc=right(i);
	// 		if (rc<(int)heap.size()) {
	// 		    if (heap[i] < heap[rc]) {
	// 			cerr<<"heap out of order"<<endl;
	// 			BREAK;
	// 		    }

	// 		    verify(rc);
	// 		}
	// 	}
	// }

	const std::vector<T>& contents(int level) const {
		return heap[level];
	}



private:



	// update the heap because the element at position i may be in the wrong place
	// usefull for pushing/popping/updating priorities
	void update_position(const int level, const int li) {

		std::vector<T> &lheap = heap[level];

		// check if we have to move up
		int p = parent(li);
		if (p>=0 && lheap[p]<lheap[li]) {
			std::swap(lheap[li], lheap[p]);
			heap_position(lheap[li], (level>0) ? -li-1 : li);
			heap_position(lheap[p],  (level>0) ? -p-1  : p);
			update_position(level, p);
		}


		// see if we have to move down
		int lc=left(li), rc=right(li);

		int maxc = -1;

		if (lc<(int)lheap.size() && rc<(int)lheap.size()) {
			maxc = (lheap[lc] < lheap[rc]) ? rc : lc;
		} else if (lc<(int)lheap.size()) {
			maxc = lc;
		} else if (rc<(int)lheap.size()) {
			maxc = rc;
		}

		// if we're smaller than the max child, swap with it
		if (maxc>=0 && lheap[li]<lheap[maxc]) {
			std::swap(lheap[li], lheap[maxc]);
			heap_position(lheap[li], (level>0) ? -li-1 : li);
			heap_position(lheap[maxc], (level>0) ? -maxc-1 : maxc);
			update_position(level, maxc);
		}
	}



	
    // zero-based heap
    static int left(int ix)    { return 2 * (ix + 1) - 1;    };
    static int right(int ix)   { return 2 * (ix + 1);        };
    static int parent(int ix)  { return ((ix + 1) / 2) - 1;  };


	std::vector<T> heap[2];

};






#endif
