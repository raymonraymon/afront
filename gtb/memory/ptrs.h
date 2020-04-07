
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


// -*-c++-*- 
/*
 * Copyright (c) 1999
 * Shachar Fleishman shacharf@math.tau.ac.il
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted 
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

/*
Auto pointers, smart pointers and auto free/close classes

  aptr - Autopointer class, A class that behaves like a pointer but frees
         it's contentens when destructed.
  aaptr - Same as aptr but for array pointers
  sptr - Smart pointer, A class that behaves like a pointer with a reference
         count. When the reference count reaches zero, the pointer if freed.
  afree - An auto free class, given an object and a function that frees
          the object, the afree class will free the object upon destruction.
  piter - An iterator adaptor that retrieves the pointer to the object.

  The aptr class sample:
  {
    aptr<C> c = new C(10);

        c->doit();
  } // c is deleted.

  Using sptr is simlilat to aptr, except that when you need another copy
  of a pointer to the same object, you use the sptr<T> instead of T*. you
  can use either the constructor or the operator= to initialize any pointer
  from the second reference to the object on.

  afree Sample:
        {
                afree<FILE*> aclose(fopen(file_name, "w"), fclose);

                ...
        } // calls fclose(f);

*/

#ifndef __PTRS_H
#define __PTRS_H

#include "counter.h"

#include <gtb/common.hpp>
#include <iterator>
#include <memory>

GTB_BEGIN_NAMESPACE

template <class T>
class aptr
{
public:
    typedef T* ptr;
    typedef const T* const_ptr;
    typedef T& reference;
    
    aptr(ptr p) : m_p(p) {}
    ~aptr() { delete m_p; }
    
    ptr operator=(ptr p) { delete m_p; return m_p = p; }
    
    bool operator==(const_ptr p) { return m_p == p; }
    
    reference operator*() { return *m_p; }
    ptr operator->() { return m_p; }
    const_ptr operator->() const { return m_p; }
    
    operator T*() { return m_p; }
    operator const T*() const { return m_p; }
    
    ptr operator+(int d) { return m_p + d; }
    ptr operator-(int d) { return m_p - d; }
    
    ptr get() { return m_p; }
    const_ptr get() const { return m_p; }
    ptr release() { ptr v = m_p; m_p = 0; return v;}
    void set(ptr p) { m_p = p; }
    
protected:
    ptr m_p;
};

template <class T>
class aaptr
{
public:
        typedef T* ptr;
        typedef const T* const_ptr;
        typedef T& reference;

        aaptr(ptr p) : m_p(p) {}
        ~aaptr() { delete[] m_p; }

        ptr operator=(ptr p) { delete[] m_p; return m_p = p; }
        bool operator==(const_ptr p) { return m_p == p; }

        reference operator*() { return *m_p; }
        ptr operator->() { return m_p; }

        operator T*() { return m_p; }
        operator const T*() const { return m_p; }

        ptr operator+(int d) { return m_p + d; }
        ptr operator-(int d) { return m_p - d; }

        ptr get() { return m_p; }
        const ptr get() const { return m_p; }
        ptr release() { ptr v = m_p; m_p = 0; return v;}
        void set(ptr p) { m_p = p; }

//        T& operator[](int p) { return m_p[p]; }
protected:
        ptr m_p;
};


template <class T>
class sptr
{
public:
    typedef T* ptr;
    typedef const T* const_ptr;
    typedef T& reference;
    
    sptr() : m_p(0), m_counter(new Counter(1)) {} // Allow storing sptr<X> in stl containers such as map
    sptr(ptr p) : m_p(p) { m_counter = new Counter(1); }
    sptr(const sptr &p) : m_p(p.m_p), m_counter(p.m_counter) { ++*m_counter; }
    ~sptr() { p_erase(); }
    
    sptr& operator=(const sptr& p) { p_erase(); m_p = p.m_p; m_counter = p.m_counter; ++*m_counter; return *this; }
    
    /*
     * Assign a new pointer, this routine creates a new sptr,
     * similar to sptr<T> p(new T);
     */
    void assign(ptr p) { p_erase(); m_counter=new Counter(1); m_p = p; }

    ptr release() 
    {
        ptr v = m_p;
        if (--*m_counter)
        {
            // BUGBUG change behaviour such that m_p==0 ==> no counter.
            m_counter = new Counter(1);
        }
        else
        {
            ++*m_counter;
        }
        m_p = 0;
        return v;
    }
    
    bool operator==(const sptr& p) const { return m_p == p.m_p; }
    bool operator==(const ptr& p) const { return m_p == p; }
    bool operator<(const sptr& p) const { return m_p < p.m_p; }
    bool operator<(const ptr& p) const { return m_p < p; }

    bool operator>(const sptr& p) const { return m_p > p.m_p; } // BUGBUG VC BUG?

    reference operator*() { return *m_p; }
    ptr operator->() { return m_p; }
    const_ptr operator->() const { return m_p; }
    
    operator T*() { return m_p; }
    operator const T*() const { return m_p; }
    
    ptr operator+(int d) { return m_p + d; }
    ptr operator-(int d) { return m_p - d; }
    
    ptr get() { return m_p; }
    const_ptr get() const { return m_p; }

protected:
    ptr m_p;
    Counter* m_counter;
    
    void p_erase() { if ( !--*m_counter) {delete m_p; delete m_counter;} }
};

template <class T>
class afree
{
public:
    typedef void (*f_free)(T);
    typedef int (*int_f_free)(T);
    typedef long (*long_f_free)(T);
    
    typedef T value_type;
    
    afree(T v, f_free free_function, T null_value = T(0) ) : 
        m_v(v), m_free(free_function) , m_null_value(null_value) {}
    afree(T v, int_f_free free_function, T null_value = T(0) ) : 
        m_v(v), m_free(reinterpret_cast<f_free>(free_function)), 
        m_null_value(null_value) {}
    afree(T v, long_f_free free_function, T null_value = T(0) ) : 
        m_v(v), m_free(reinterpret_cast<f_free>(free_function)),
        m_null_value(null_value)  {}
    ~afree() { if (m_v != m_null_value) m_free(m_v); }
    
    operator T&() { return m_v; }
    T& operator ->() { return m_v; }
protected:
    T m_v;
    f_free m_free;
    T m_null_value;
};

//
// TEMPLATE CLASS auto_ptr
// Taken from VC STL, release keeps the pointer.
//
// Allows me to store pointers in containers
// And... mix owned and non-owned pointers in the one
// container.
//
template<class T>
class auto_ptr2 {
public:
    typedef T element_type;

    explicit auto_ptr2(element_type *p = 0) throw() :
        m_owner(p != 0), m_ptr(p)
    {
    }

    auto_ptr2(const auto_ptr2& p) throw():
        m_owner(p.m_owner), m_ptr(p.release())
    {
    }    

    auto_ptr2& operator=(const auto_ptr2& rhs) throw()
    {
        if (this != &rhs)
        {
            if (m_ptr != rhs.get())
            {
                if (m_owner) delete m_ptr;
                m_owner = rhs.m_owner;
            }
            else if (rhs.m_owner) m_owner = true;
            m_ptr = rhs.release();
        }
        return *this; 
    }

    ~auto_ptr2()
    {
        if (m_owner) delete m_ptr;
    }

    element_type& operator*() const throw()
    {
        return *get();
    }

    element_type *operator->() const throw()
    {
        return get();
    }

    element_type *get() const throw()
    {
        return m_ptr;
    }

    element_type *release() const throw()
    {
        m_owner = false;
        return m_ptr;
    }
private:
    mutable bool m_owner;
    element_type *m_ptr;
}; //auto_ptr2

//
// Iterator adaptor that retrieves the pointer to the object
//
template <class CONT>
class piter : public std::iterator<std::random_access_iterator_tag,
				   typename CONT::value_type,
				   int>
{
public:
    typedef typename CONT::value_type* value_type;
    typedef typename CONT::iterator iterator;
    iterator i;

    piter() {}
    piter(iterator iter) : i(iter) {}

    typename CONT::value_type* operator*() { return &(*i); }
    piter& operator++() { ++i; return *this; }
    int operator-(const piter& rhs) const { return i - rhs.i; }

    bool operator ==(const piter& rhs) { return i == rhs.i; }
    bool operator !=(const piter& rhs) { return i != rhs.i; }

    piter& operator=(const piter& rhs) { i = rhs.i; return *this; }
    piter& operator=(iterator rhs) { i = rhs;  return *this; }

};
GTB_END_NAMESPACE

#endif // __PTRS_H
