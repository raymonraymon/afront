
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


#ifndef __KNN_H
#define __KNN_H

/*
 * A data structure to find fast the K nn of a set of points
 * The data structure is a 2D array of NxK of point indices
 * which point to surfelset
 *
 * Support:
 *    - Creation
 *    - Initializing the data strucutre by inserting the neighbors
 *      of a point
 *    - query
 *
 * Resizing is not supported!!!
 *
 * N - number of points
 * K - number of nearest neighbors / point
 *
 * Support routines
 *   Generate from pointset
 *
 * Note: Not thread-safe
 */

#include "mlslibdefs.h"
#include <gtb/gtb.hpp>

//#include <threadslib/threadslib.h>

MLSLIB_BEGIN_NAMESPACE

template <class REAL>
class KNN
{
public:
    typedef unsigned index_type;
    typedef gtb::tsurfel_set<REAL> surfel_set;
    typedef gtb::tsurfelset_view<REAL> surfelset_view;
    typedef gtb::tPoint3<REAL> Point3;

    KNN(unsigned K, surfel_set& surfels);
    virtual ~KNN();

    void set(index_type idx, index_type k, index_type value);
    index_type get(index_type idx, index_type k) const;

    void extract_connectivity(index_type origin, unsigned radius, surfelset_view& view);
    void extract(index_type origin, REAL radius, surfelset_view& view);
    void extract(index_type origin, unsigned K, surfelset_view& view);

    surfel_set& get_points();

    unsigned getK() const { return _K; }
    unsigned getN() const { return _N; }

protected:
    unsigned _N;
    unsigned _K;

    gtb::aaptr<index_type> _narray;


    // The following is used do avoid extracting a point twice
//    unsigned last_extract_value;
//    std::vector<unsigned> t_visited_iteration;
    typedef std::vector<unsigned> t_visited_iteration;

    //
    // Thead-safe last_extract value and visited_iteration
    // [ Windows specific ]
    //
    typedef std::vector<unsigned> t_tls_last_extract_value;
    typedef std::list<t_visited_iteration> t_tls_visited_iteration;

    t_tls_last_extract_value _tls_last_extract_value;
    t_tls_visited_iteration _tls_visited_iteration;

    unsigned increase_last_extract_value();
    t_visited_iteration& get_visited_iteration();
    unsigned get_my_thread_index();
//    thlib::CSObject _guard;
    unsigned _num_threads;

    long _tls_index;
    // [ Thead-safe ]

    struct extract_node
    {
        extract_node() {}
        extract_node(index_type index, unsigned radius) : _index(index), _radius(radius) {}

        index_type _index;
        unsigned _radius;
    };

    typedef std::queue<extract_node> t_stack;

    surfel_set& _surfels;
};

template <class REAL>
inline void KNN<REAL>::set(index_type idx, index_type k, index_type value)
{
    assert(idx < _N);
    assert(k < _K);
    _narray[idx*_K+k] = value;
}

template <class REAL>
inline typename KNN<REAL>::index_type KNN<REAL>::get(index_type idx, index_type k) const
{
    assert(idx < _N);
    assert(k < _K);

    return _narray[idx*_K+k];
}

template <class REAL>
void build_knn(KNN<REAL>& knn);

MLSLIB_END_NAMESPACE

#endif // __KNN_H
