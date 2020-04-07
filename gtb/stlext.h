
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


#ifndef __STLEXT_H
#define __STLEXT_H


// compact a vector since erase never frees any memory
template <typename T>
void compact(std::vector<T> &v) {
	if (v.size() != v.capacity()) {
		typedef std::vector<T> vecT;
		vecT tmp = v;

		// clear v's memory
		v.~vecT();
		new (&v) vecT;

		v = tmp;
	}
}




// Need to find the right place for this as well


/*
 * AssignOp - a function object that assigns a value to an iterator
 *      sample use in kdtree.h:         Traverse(point, K, radius, AssignOp<T>(it));
 *
 * mem_fun2 - same as mem_fun1 but for member functions that receive
 *            two parameters
 *
 * counter_iterator - an iterator that generates a sequence of number: 1,2,3...
 *                    for example set_diff(counter_iterator(0), counter_iterator(10), i1, i2, dest);
 *
 * func_composite - an adaptor that gets two functions f,g and for each argument x returns f(g(x))
 */

#include <functional>

#if 0
template<class _Bfn>
class void_binder1st : 
public std::unary_function<typename _Bfn::second_argument_type, typename _Bfn::result_type> 
{
public:
	void_binder1st(const _Bfn& _X,
		const typename _Bfn::first_argument_type& _Y)
		: op(_X), value(_Y) {}
	void operator()(const typename _Bfn::second_argument_type& _X) const
	{
        op(value, _X); 
    }
protected:
	_Bfn op;
	_Bfn::first_argument_type value;
	};
#else
template <class _Operation>
class void_binder1st
  : public std::unary_function<typename _Operation::second_argument_type, void> {
protected:
  _Operation op;
  typename _Operation::first_argument_type value;
public:
  void_binder1st(const _Operation& __x,
            const typename _Operation::first_argument_type& __y)
	: op(__x), value(__y) {}
  void
  operator()(const typename _Operation::second_argument_type& __x) const {
    op(value, __x);
  }
};

#endif

template<class _Bfn, class _Ty> inline
void_binder1st<_Bfn> void_bind1st(const _Bfn& XX, const _Ty& YY)
{
    return (void_binder1st<_Bfn>(XX, _Bfn::first_argument_type(YY))); 
}

template<class _Ty, class AA>
class void_mem_fun1_t : public std::binary_function<_Ty*, AA, void>
{
private:
	void (_Ty::*_Ptr)(AA);

public:
	explicit void_mem_fun1_t(void (_Ty::*_Pm)(AA)) : _Ptr(_Pm) {}
	void operator()(_Ty *PP, AA _Arg) const
    {
        (PP->*_Ptr)(_Arg); 
    }
};

template<class _Ty, class AA> 
inline void_mem_fun1_t<_Ty, AA> void_mem_fun1(void (_Ty::*_Pm)(AA))
{
    return (void_mem_fun1_t<_Ty, AA>(_Pm)); 
}

//
// A test: memfun that works on objects and not pointers to objects
//
// The problem
//   if you have something like:
//     vector<myclass> mc;
//   you cannot
//     for_each(mc.begin(), mc.end(), mem_fun(&mc::doit));
//   because:
//     mem_fun_t::operator() (mc* x) and not mem_fun_t::operator() (mc& x) !!!
//
		// TEMPLATE CLASS const_mem_fun_t
template<class _Result,
	class _Ty>
	class const_mem_fun_t_np
        : public std::unary_function<const _Ty *, _Result>
	{	// functor adapter (*p->*pfunc)(), const *pfunc
public:
	explicit const_mem_fun_t_np(_Result (_Ty::*_Pm)() const)
		: _Pmemfun(_Pm)
		{	// construct from pointer
		}

	_Result operator()(const _Ty& _Pleft) const
		{	// call function
		return (((&_Pleft)->*_Pmemfun)());
		}

private:
	_Result (_Ty::*_Pmemfun)() const;	// the member function pointer
	};

template<class _Result,
	class _Ty> inline
	const_mem_fun_t_np<_Result, _Ty>
		mem_fun_np(_Result (_Ty::*_Pm)() const)
	{	// return a const_mem_fun_t functor adapter
	return (const_mem_fun_t_np<_Result, _Ty>(_Pm));
	}

/*---------------------------------*/
/*
 *  mem_fun2:
 *    for a function to call a member that received two arguments
 */

template<class _Result,
         class _Ty, // Class Type
         class _Arg1, // first argument type
         class _Arg2> // Second argument type
class mem_fun2_t
{	// functor adapter (*p->*pfunc)(val1, val2), non-const *pfunc
public:
    explicit mem_fun2_t(
        _Ty* object,
        _Result (_Ty::*_Pm)(_Arg1, _Arg2))
        : _object(object), _Pmemfun(_Pm)
    {	// construct from pointer
    }

    _Result operator()(_Arg1 a1, _Arg2 a2) const
    {	// call function with operand
        return ((_object->*_Pmemfun)(a1, a2));
    }

private:
    _Ty* _object;
    _Result (_Ty::*_Pmemfun)(_Arg1, _Arg2);	// the member function pointer
};

/*
 * the ""mem_fun2"" function
 */
template<class _Result,
         class _Ty,
         class _Arg1,
         class _Arg2>
inline mem_fun2_t<_Result, _Ty, _Arg1,_Arg2> mem_fun2(_Ty* object, _Result (_Ty::*_Pm)(_Arg1, _Arg2))
{	// return a mem_fun2_t functor adapter
    return (mem_fun2_t<_Result, _Ty, _Arg1, _Arg2>(object, _Pm));
}


template<class _Result,
         class _Ty, // Class Type
         class _Arg1, // first argument type
         class _Arg2> // Second argument type
class const_mem_fun2_t
{	// functor adapter (*p->*pfunc)(val1, val2), non-const *pfunc
public:
    explicit const_mem_fun2_t(
        _Ty* object,
        _Result (_Ty::*_Pm)(_Arg1, _Arg2) const)
        : _object(object), _Pmemfun(_Pm)
    {	// construct from pointer
    }

    _Result operator()(_Arg1 a1, _Arg2 a2) const
    {	// call function with operand
        return ((_object->*_Pmemfun)(a1, a2));
    }

private:
    _Ty* _object;
    _Result (_Ty::*_Pmemfun)(_Arg1, _Arg2) const;	// the member function pointer
};

template<class _Result,
         class _Ty, // Class Type
         class _Arg1, // first argument type
         class _Arg2,// Second argument type
         class _Arg3>  // THIRD AT
class const_mem_fun3_t
{	// functor adapter (*p->*pfunc)(val1, val2), non-const *pfunc
public:
    explicit const_mem_fun3_t(
        _Ty* object,
        _Result (_Ty::*_Pm)(_Arg1, _Arg2, _Arg3) const)
        : _object(object), _Pmemfun(_Pm)
    {	// construct from pointer
    }

    _Result operator()(_Arg1 a1, _Arg2 a2, _Arg3 a3) const
    {	// call function with operand
        return ((_object->*_Pmemfun)(a1, a2, a3));
    }

private:
    _Ty* _object;
    _Result (_Ty::*_Pmemfun)(_Arg1, _Arg2, _Arg3) const;	// the member function pointer
};

/*--------------- function object AssignOP ----------------*/
//
// Assign value to an iterator

template<class ITERATOR, class T>
struct AssignOP_t : public std::unary_function<T, void>
{
    mutable ITERATOR it;

    AssignOP_t(ITERATOR IT) : it(IT) {}

    void operator() (const T& x) const
    {
        *it = x;
    }
};

template<class T, class ITERATOR>
AssignOP_t<ITERATOR, T> AssignOp(ITERATOR it)
{
    return AssignOP_t<ITERATOR, T>(it);
}

/*--------------- counter_iterator -------------------*/
class counter_iterator : public std::random_access_iterator_tag
{
public:
    counter_iterator(const counter_iterator& rhs) : _v(rhs._v) {}
    counter_iterator(int value) : _v(value) {}

    bool operator==(const counter_iterator& rhs) const { return _v == rhs._v; }
    bool operator!=(const counter_iterator& rhs) const { return !(*this == rhs); }
    bool operator<(const counter_iterator& rhs) const { return _v < rhs._v; }

    counter_iterator& operator++() { ++_v; return *this; }
    counter_iterator operator++(int) {counter_iterator r(*this); ++_v; return r; }
    counter_iterator& operator+=(const counter_iterator& rhs) { _v += rhs._v; return *this; }
    counter_iterator& operator+=(int delta) { _v += delta; return *this; }
    counter_iterator operator+(const counter_iterator& rhs) { counter_iterator r(*this); r += rhs; return r; }
    counter_iterator operator+(int delta) { counter_iterator r(*this); r += delta; return r; }

    counter_iterator& operator--() { --_v; return *this; }
    counter_iterator operator--(int) {counter_iterator r(*this); --_v; return r; }
    counter_iterator& operator-=(const counter_iterator& rhs) { _v -= rhs._v; return *this; }
    counter_iterator& operator-=(int delta) { _v -= delta; return *this; }
    counter_iterator operator-(const counter_iterator& rhs) { counter_iterator r(*this); r -= rhs; return r; }
    counter_iterator operator-(int delta) { counter_iterator r(*this); r -= delta; return r; }

    int operator*() const { return _v; }

protected:
    int _v;
};


/*----------------- std::pair compare function objects -----------------*/
/*
 * Function objects to compare the first or second members of std::pari
 */
template<class T1, class T2>
struct s_stdpaircompare1
{
    bool operator () (const std::pair<T1,T2>& lhs, const std::pair<T1,T2>& rhs) const { return lhs.first < rhs.first; }
};

template<class T1, class T2>
struct s_stdpaircompare2
{
    bool operator () (const std::pair<T1,T2>& lhs, const std::pair<T1,T2>& rhs) const { return lhs.second < rhs.second; }
};

template<class T>
s_stdpaircompare1<typename T::value_type::first_type, typename T::value_type::second_type> stdpaircompare1(const T& cont)
{
    s_stdpaircompare1<typename T::value_type::first_type, typename T::value_type::second_type> cmp;
    return cmp;
}

template<class T>
s_stdpaircompare2<typename T::value_type::first_type, typename T::value_type::second_type> stdpaircompare2(const T& cont)
{
    s_stdpaircompare2<typename T::value_type::first_type, typename T::value_type::second_type> cmp;
    return cmp;
}

#endif // __STLEXT_H
