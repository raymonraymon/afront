
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


#include "stdafx.h"
#include "knn.h"

MLSLIB_BEGIN_NAMESPACE

template <class REAL>
KNN<REAL>::KNN(unsigned K, surfel_set& surfels) :
    _N(surfels.size()),
    _K(K),
    _narray(new index_type[_N*K]),
    _num_threads(0),
    _surfels(surfels)
{
  /*
    _tls_index = TlsAlloc();
    if (_tls_index == TLS_OUT_OF_INDEXES)
    {
        printf("KNN::KNN failed to allocate tls\n");
    }
  */
}

template <class REAL>
KNN<REAL>::~KNN()
{
  /*
    TlsFree(_tls_index);
  */
}

/*
 * Extract nearest neighbors to some "origin"
 * using the connectivity, going "radius" neighborhood steps in the 
 * connectivity map
 */
template <class REAL>
void KNN<REAL>::extract_connectivity(index_type origin, unsigned radius, surfelset_view& view)
{
    unsigned last_extract_value = increase_last_extract_value();
    t_visited_iteration& _visited_iteration = get_visited_iteration();

    _visited_iteration[origin] = last_extract_value;
    view.insert(origin);

    t_stack s;
    s.push(extract_node(origin, radius));
    while (!s.empty())
    {
        extract_node node = s.front();
        s.pop();

        index_type idx = node._index;
        unsigned next_radius = node._radius-1;

        for (unsigned i = 0 ; i < _K; ++i)
        {
            index_type neighbor_idx = _narray[idx*_K+i];
            if (_visited_iteration[neighbor_idx] != last_extract_value)
            {
                _visited_iteration[neighbor_idx] = last_extract_value;
                view.insert(neighbor_idx);
                if (next_radius) s.push(extract_node(neighbor_idx, next_radius));
            }
        }
    }
}

template <class REAL>
void KNN<REAL>::extract(index_type origin, REAL radius, surfelset_view& view)
{
    unsigned last_extract_value = increase_last_extract_value();
    t_visited_iteration& _visited_iteration = get_visited_iteration();

    _visited_iteration[origin] = last_extract_value;
    view.insert(origin);

    REAL radius2 = radius * radius;

    const Point3& origin3 = _surfels.vertex(origin);

    std::queue<unsigned> s;
    s.push(origin);

    while (!s.empty())
    {
        unsigned idx = s.front();
        s.pop();

        // const Point3& p = _surfels.vertex(idx); Unused

        for (unsigned i = 0 ; i < _K; ++i)
        {
            index_type neighbor_idx = _narray[idx*_K+i];
            if (_visited_iteration[neighbor_idx] != last_extract_value)
            {
                _visited_iteration[neighbor_idx] = last_extract_value;

                REAL d2 = (_surfels.vertex(neighbor_idx) - origin3).squared_length();

                if (d2 <= radius2)
                {
                    view.insert(neighbor_idx);
                    s.push(neighbor_idx);
                }
            }
        }
    }
}

//
// Extract K nearest neighbors
//
template <class REAL>
void KNN<REAL>::extract(index_type origin, unsigned K, surfelset_view& view)
{
    unsigned last_extract_value = increase_last_extract_value();
    t_visited_iteration& _visited_iteration = get_visited_iteration();

    _visited_iteration[origin] = last_extract_value;

    const Point3& origin3 = _surfels.vertex(origin);

    //
    // a map from distance to point index
    //
    std::multimap<REAL, unsigned> s;
    typedef typename std::multimap<REAL, unsigned>::value_type mapval;

    s.insert(mapval(0.0, origin));

    while (!s.empty())
    {
        unsigned idx = s.begin()->second;
        s.erase(s.begin());
        view.insert(idx);
        if (view.size() >= K) return;

        // const Point3& p = _surfels.vertex(idx); Unused variable

        for (unsigned i = 0 ; i < _K; ++i)
        {
            index_type neighbor_idx = _narray[idx*_K+i];
            if (_visited_iteration[neighbor_idx] != last_extract_value)
            {
                _visited_iteration[neighbor_idx] = last_extract_value;

                REAL d2 = (_surfels.vertex(neighbor_idx) - origin3).squared_length();
                s.insert(mapval(d2, neighbor_idx));
            }
        }
    }
}

template <class REAL>
typename KNN<REAL>::surfel_set& KNN<REAL>::get_points()
{
    return _surfels;
}

/*------ threads stuff --------*/
template <class REAL>
unsigned KNN<REAL>::increase_last_extract_value()
{
//    thlib::CS __cs(_guard);
    unsigned my_index = get_my_thread_index();
    ++(_tls_last_extract_value[my_index]);
    return _tls_last_extract_value[my_index];
}

template <class REAL>
typename KNN<REAL>::t_visited_iteration& KNN<REAL>::get_visited_iteration()
{
//    thlib::CS __cs(_guard);
    unsigned my_index = get_my_thread_index();
    t_tls_visited_iteration::iterator p_visited = _tls_visited_iteration.begin();
    std::advance(p_visited, my_index);
    return *p_visited;
}

/*
 * Return my thread index
 * if the current thread was never used before, allocate an index to it
 * and initialize its last_value and visited_iteration
 */
template <class REAL>
unsigned KNN<REAL>::get_my_thread_index()
{
	return 0;
	/*
    unsigned my_index = (unsigned)TlsGetValue(_tls_index);
    if (my_index == 0)
    {
//        thlib::CS __cs(_guard);
        ++_num_threads;
        my_index = _num_threads; // store 1 for the first thread
        TlsSetValue(_tls_index, (LPVOID)my_index);
        
        _tls_last_extract_value.push_back(0);
        t_visited_iteration varray;
        _tls_visited_iteration.push_back(varray);
        _tls_visited_iteration.back().resize(_N, 0);

        return my_index-1;
    }
    else
    {
        return my_index-1; // we want to start from 0, but if TlsGetValue returns
                           // a 0, it means an error, so we start with one
    }
	*/
}

/*------ [ threads stuff ] --------*/


// Instantiate KNN
template class KNN<float>;
template class KNN<double>;

/*------------- Tools etc -----------------*/


template <class REAL>
void build_knn(KNN<REAL>& knn)
{
    typedef gtb::tsurfel_set<REAL> surfel_set;
    typedef gtb::tsurfelset_view<REAL> surfelset_view;
    typedef gtb::KDTree<int, REAL, gtb::GetPoint_f<surfel_set> > t_surfel_tree;

    surfel_set& points = knn.get_points();
    unsigned N = points.size();
    unsigned K = knn.getK();

    // Generate the tree
    gtb::aptr<t_surfel_tree> tree = new t_surfel_tree(5, points.bounding_box(), gen_GetPoint_f(points));
    for (unsigned i = 0; i < N; ++i)
    {
        tree->Insert(i);
    }
    tree->MakeTree();

    // Generate the KNN, by query
    surfelset_view neighbors(points);
    for (unsigned i = 0; i < N; ++i)
    {
        neighbors.clear();
        const gtb::tPoint3<REAL>& p = points.vertex(i);
        tree->Extract(p, K+1, 1e8, std::back_inserter(neighbors.get_view()));
        for (unsigned j = 0; j < K; ++j) // BUG? 1..K and not 0..K-1??
        {
            knn.set(i, j, neighbors.get_index(j+1));
        }

        if ((i % 100000) == 0)
        {
            printf("KNN: %d%%    \r", i*100/N);
        }
    }
    printf("                                                                \r");
}

template void build_knn(KNN<float>& knn);
template void build_knn(KNN<double>& knn);


MLSLIB_END_NAMESPACE
