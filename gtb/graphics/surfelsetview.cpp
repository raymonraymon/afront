
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
#include <gtb/graphics/surfelset.hpp>
#include <gtb/graphics/surfelsetview.hpp>

GTB_BEGIN_NAMESPACE

template <class T>
tsurfelset_view<T>::tsurfelset_view(const tsurfelset_view& rhs) :
    _surfelset(rhs._surfelset),
    _view(rhs._view)
{
}

template <class T>
tsurfelset_view<T>::tsurfelset_view(const tsurfel_set<T>& surfelset) :
    _surfelset(surfelset)
{
}


/*
 * Find the position where point with given index.
 * BAD function: linear time
 */
template <class T>
unsigned tsurfelset_view<T>::find_index(unsigned idx) const
{
	for (unsigned i = 0; i < size(); ++i)
	{
		if (_view[i] == idx) return i;
	}
	
#ifndef NO_EXCEPTIONS
	throw CErr("tsurfelset_view::find_index: not found");
#else
	return -1;
#endif
}

template <class T>
void tsurfelset_view<T>::insertall()
{
	int N = _surfelset.size();
	for (int i = 0; i < N; ++i)
	{
		_view.push_back(i);
	}
}

template <class T>
void tsurfelset_view<T>::compute_bounding_box() const
{
    int N = size();
    tPoint3<T> r(-1e8,-1e8,-1e8);
    tPoint3<T> l(1e8,1e8,1e8);
    for (int i = 0; i < N; ++i)
    {
        const tPoint3<T>& xi = vertex(i);
        for (int j = 0; j < 3; ++j)
        {
            l[j] = min2(l[j], xi[j]);
            r[j] = max2(r[j], xi[j]);
        }
    }
    tModel<T>::_bounding_box = tBox3<T>(l,r);
}

template <class T>
void tsurfelset_view<T>::compute_centroid() const
{
    tModel<T>::_centroid = tPoint3<T>(0,0,0);
    int N = size();
    for (int i = 0; i < N; ++i)
    {
        const tPoint3<T>& xi = vertex(i);
        for (int j = 0; j < 3; ++j)
        {
            tModel<T>::_centroid[j] += xi[j];
        }
    }
    if (N>0) tModel<T>::_centroid.scalar_scale(1.0 / N);
}

// the median of the set of points in each axis indipendently...
template <class T>
void tsurfelset_view<T>::compute_median() const
{
    int N = size();
    std::vector<T> values;
    values.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        values.push_back(vertex(i).x());
    }
    std::sort(values.begin(), values.end());
    tModel<T>::_median[0] = values[N/2];

    for (int i = 0; i < N; ++i)
    {
        values[i] = vertex(i).y();
    }
    std::sort(values.begin(), values.end());
    tModel<T>::_median[1] = values[N/2];

    for (int i = 0; i < N; ++i)
    {
        values[i] = vertex(i).z();
    }
    std::sort(values.begin(), values.end());
    tModel<T>::_median[2] = values[N/2];
}

template <class T>
void tsurfelset_view<T>::fill(int start, int end)
{
    assert(start <= end);
    for (int i = start; i <= end; ++i)
    {
        insert(i);
    }
}

template <class T>
tsurfelset_view<T>& tsurfelset_view<T>::operator=(const tsurfelset_view& rhs)
{
    assert(&_surfelset == &rhs._surfelset);
    clear();
    std::copy(rhs._view.begin(), rhs._view.end(), std::back_inserter(_view));
    return *this;
}

/*
 * Erase the indices from this view that are in the rhs view
 * Side effect: The order of the view changes!!!!
 */
template <class T>
void tsurfelset_view<T>::erase(tsurfelset_view& rhs)
{
    std::sort(rhs._view.begin(), rhs._view.end());

    subset_view sorted_view(_view);
    std::sort(sorted_view.begin(), sorted_view.end());
    subset_view new_view;
    std::set_difference(
        sorted_view.begin(), 
        sorted_view.end(),
        rhs._view.begin(),
        rhs._view.end(),
        std::back_inserter(new_view));
    _view = new_view;
}

template <class T>
void tsurfelset_view<T>::erase(unsigned idx)
{
    assert(idx < _view.size());
    _view.erase(_view.begin()+idx);
}

template <class T>
void tsurfelset_view<T>::erase(unsigned idx1, unsigned idx2)
{
    assert(idx1 < _view.size());
    assert(idx1 < _view.size());
    _view.erase(_view.begin()+idx1, _view.begin()+idx2);
}

/*
 * Print the list of indices
 */
template <class T>
void tsurfelset_view<T>::printme(const char* header)
{
    printf("%s", header);
    for (unsigned jj = 0; jj < size(); ++jj)
    {
        printf(" %d", get_index(jj));
    }
    printf("\n");
}

template class tsurfelset_view<float>;
template class tsurfelset_view<double>;

GTB_END_NAMESPACE
