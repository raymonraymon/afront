
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


#ifndef __SEQUENCEVIEWCLASS_H
#define __SEQUENCEVIEWCLASS_H

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

/*
 * A general purpose view class to an array of object
 * A view to an array is behaves like an array to a subset of
 * a "main" array. (Is that clear??)
 *
 * In other words, this class holds a list of indices to the main array
 * and allows referencing only the subset of objects in the main
 * array that are in the list of indices
 *
 *
 * The template argument is for the type of the parent (main) array
 *
 * For laziness reasons, the parent class is not const!!!
 *
 * NEW: support sorting the view with a function object that compares
 *      members of the container items
 */

template <class T>
class SequenceView
{
public:
    typedef typename T::value_type value_type;
    typedef typename T::reference reference;
    typedef typename T::const_reference const_reference;

    typedef T parent_type;
    typedef std::vector<unsigned> indices_type;

    class const_iterator : public std::iterator<std::random_access_iterator_tag, const value_type>
    {
    public:
        typedef typename T::reference reference;
        typedef typename T::pointer pointer;
	typedef typename T::difference_type difference_type;

        const_iterator();
        const_iterator(const const_iterator& rhs);
        const_iterator(indices_type::iterator idx, parent_type* parent);

        reference operator*() const;
        pointer operator->() const;
        const_iterator& operator++();
        const_iterator operator++(int);
        const_iterator& operator--();
        const_iterator operator--(int);
        const_iterator& operator+=(difference_type _Off);
        const_iterator operator+(difference_type _Off) const;
        const_iterator& operator-=(difference_type _Off);
        const_iterator operator-(difference_type _Off) const;
        difference_type operator-(const const_iterator& _Right) const;
        reference operator[](difference_type _Off) const;

        bool operator==(const const_iterator& _Right) const;
        bool operator!=(const const_iterator& _Right) const;
        bool operator<(const const_iterator& _Right) const;
        bool operator>(const const_iterator& _Right) const;
        bool operator<=(const const_iterator& _Right) const;
        bool operator>=(const const_iterator& _Right) const;

    protected:
        parent_type* _parent;
        indices_type::iterator _idx;
        friend class SequenceView<T>;
    };

    class iterator : public const_iterator
    {
    public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef typename T::value_type value_type;
        typedef typename T::reference reference;
        typedef typename T::pointer pointer;
	typedef typename T::difference_type difference_type;

        iterator();
        iterator(const iterator& rhs);
        iterator(indices_type::iterator idx, parent_type* parent);

        reference operator*() const;
        pointer operator->() const;
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        iterator& operator+=(difference_type _Off);
        iterator operator+(difference_type _Off) const;
        iterator& operator-=(difference_type _Off);
        iterator operator-(difference_type _Off) const;
        difference_type operator-(const iterator& _Right) const;
        reference operator[](difference_type _Off) const;

    protected:
        parent_type* _parent;
        indices_type::iterator _idx;
        friend class SequenceView<T>;
    };

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    typedef indices_type::size_type size_type;

    SequenceView(parent_type& parent);

// vector like behaviour
    reference operator[](size_type idx);
    const_reference operator[](size_type idx) const;

    reference at(size_type idx);
    const_reference at(size_type idx) const;
    reference back( );
    const_reference back( ) const;
    const_iterator begin( ) const;
    iterator begin( );
    void clear( );
    bool empty( ) const;
    iterator end( );
    const_iterator end( ) const;
    iterator erase(iterator where);
    iterator erase(iterator first, iterator last);
    reference front( );
    const_reference front( ) const;

    reverse_iterator rbegin( );
    const_reverse_iterator rbegin( ) const;
    const_reverse_iterator rend( ) const;
    reverse_iterator rend( );

    size_type size() const;

//    void push_back(indices_type::value_type idx) { insert_index(idx); }

    // special stuff
    void insert_index(indices_type::value_type idx);
    indices_type& get_indices_vector( );
    parent_type& get_parent( );
    unsigned index(int idx); // return the idx'th index

    // make the view a see the entire parent
    void use_all();

protected:
    parent_type& _parent;
    indices_type _indices;

    //
    // A function object that is used to sort
    // the view.
    //
    template<class COMP>
    struct CompareObject
    {
        CompareObject(parent_type& container, COMP& compare) :
            _container(container), _compare(compare)
        {
        }

        bool operator()(unsigned lhs, unsigned rhs)
        {
            return _compare(_container[lhs], _container[rhs]);
        }

        parent_type _container;
        COMP& _compare;
    };

public:

    // JS FIXME
    template<class COMP>
      //    void sort(COMP& compare_func = std::less<typename parent_type::value_type>)
    void sort(COMP& compare_func)
    {
      CompareObject<COMP> cobj(_parent, compare_func);
      std::sort(_indices.begin(), _indices.end(), cobj);
    }
};

GTB_END_NAMESPACE

#include "view.inl"

#endif // __SEQUENCEVIEWCLASS_H
