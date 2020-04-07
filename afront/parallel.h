
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



#ifndef AFRONT_PARALLEL_H
#define AFRONT_PARALLEL_H


/*
some examples


-------------
// to call a member function in parallel

class someclass {
  void foo() {
    ParallelExecutor(idealNumThreads, makeClassFunctor(this, &someclass::bar));
  }
  void bar(int nt, int id) {
    // do some parallel work
  }
};

-------------
// to call a global function / functor

void pe2(int nt, int id, int bi, const char *s) {
  // do some parallel work
}

int bi=0;
char *s="";
ParallelExecutor(idealNumThreads, pe2, bi, s);



// To add more versions for more parameters to your functions, you have to define the appropriate macro's below,
// and instantiate them at the bottom.  For the class functors to work, you need 2 extra variables since the 
// id and numThreads are not treated as special.  That's why the #defines go up to 7, but the parallel functions
// can only take 5 variables.


*/



#define PE_TEMPLATE_DEF0
#define PE_TEMPLATE_DEF1 <typename T1>
#define PE_TEMPLATE_DEF2 <typename T1, typename T2>
#define PE_TEMPLATE_DEF3 <typename T1, typename T2, typename T3>
#define PE_TEMPLATE_DEF4 <typename T1, typename T2, typename T3, typename T4>
#define PE_TEMPLATE_DEF5 <typename T1, typename T2, typename T3, typename T4, typename T5>
#define PE_TEMPLATE_DEF6 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
#define PE_TEMPLATE_DEF7 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
#define PE_TEMPLATE_DEF8 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
#define PE_TEMPLATE_DEF9 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
#define PE_TEMPLATE_DEF10 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
#define PE_TEMPLATE_DEF11 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>

#define PE_TEMPLATE_ADD_DEF0
#define PE_TEMPLATE_ADD_DEF1 , typename T1
#define PE_TEMPLATE_ADD_DEF2 , typename T1, typename T2
#define PE_TEMPLATE_ADD_DEF3 , typename T1, typename T2, typename T3
#define PE_TEMPLATE_ADD_DEF4 , typename T1, typename T2, typename T3, typename T4
#define PE_TEMPLATE_ADD_DEF5 , typename T1, typename T2, typename T3, typename T4, typename T5
#define PE_TEMPLATE_ADD_DEF6 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6
#define PE_TEMPLATE_ADD_DEF7 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7
#define PE_TEMPLATE_ADD_DEF8 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8
#define PE_TEMPLATE_ADD_DEF9 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9
#define PE_TEMPLATE_ADD_DEF10 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10
#define PE_TEMPLATE_ADD_DEF11 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11

#define PE_TEMPLATE_ADD_CALL0 
#define PE_TEMPLATE_ADD_CALL1 , T1
#define PE_TEMPLATE_ADD_CALL2 , T1, T2
#define PE_TEMPLATE_ADD_CALL3 , T1, T2, T3
#define PE_TEMPLATE_ADD_CALL4 , T1, T2, T3, T4
#define PE_TEMPLATE_ADD_CALL5 , T1, T2, T3, T4, T5
#define PE_TEMPLATE_ADD_CALL6 , T1, T2, T3, T4, T5, T6
#define PE_TEMPLATE_ADD_CALL7 , T1, T2, T3, T4, T5, T6, T7
#define PE_TEMPLATE_ADD_CALL8 , T1, T2, T3, T4, T5, T6, T7, T8
#define PE_TEMPLATE_ADD_CALL9 , T1, T2, T3, T4, T5, T6, T7, T8, T9
#define PE_TEMPLATE_ADD_CALL10 , T1, T2, T3, T4, T5, T6, T7, T8, T9, T10
#define PE_TEMPLATE_ADD_CALL11 , T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11

#define PE_TEMPLATE_CALL0
#define PE_TEMPLATE_CALL1 <T1>
#define PE_TEMPLATE_CALL2 <T1,T2>
#define PE_TEMPLATE_CALL3 <T1,T2,T3>
#define PE_TEMPLATE_CALL4 <T1,T2,T3,T4>
#define PE_TEMPLATE_CALL5 <T1,T2,T3,T4,T5>
#define PE_TEMPLATE_CALL6 <T1,T2,T3,T4,T5,T6>
#define PE_TEMPLATE_CALL7 <T1,T2,T3,T4,T5,T6,T7>
#define PE_TEMPLATE_CALL8 <T1,T2,T3,T4,T5,T6,T7,T8>
#define PE_TEMPLATE_CALL9 <T1,T2,T3,T4,T5,T6,T7,T8,T9>
#define PE_TEMPLATE_CALL10 <T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>
#define PE_TEMPLATE_CALL11 <T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>

#define PE_PARAMETER_ADD_DEF0(T,p)
#define PE_PARAMETER_ADD_DEF1(T,p) ,T##1 p##1
#define PE_PARAMETER_ADD_DEF2(T,p) ,T##1 p##1, T##2 p##2
#define PE_PARAMETER_ADD_DEF3(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3
#define PE_PARAMETER_ADD_DEF4(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4
#define PE_PARAMETER_ADD_DEF5(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5
#define PE_PARAMETER_ADD_DEF6(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6
#define PE_PARAMETER_ADD_DEF7(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7
#define PE_PARAMETER_ADD_DEF8(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8
#define PE_PARAMETER_ADD_DEF9(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9
#define PE_PARAMETER_ADD_DEF10(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10
#define PE_PARAMETER_ADD_DEF11(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11

#define PE_PARAMETER_CALL0(p)
#define PE_PARAMETER_CALL1(p) p##1
#define PE_PARAMETER_CALL2(p) p##1, p##2
#define PE_PARAMETER_CALL3(p) p##1, p##2, p##3
#define PE_PARAMETER_CALL4(p) p##1, p##2, p##3, p##4
#define PE_PARAMETER_CALL5(p) p##1, p##2, p##3, p##4, p##5
#define PE_PARAMETER_CALL6(p) p##1, p##2, p##3, p##4, p##5, p##6
#define PE_PARAMETER_CALL7(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7
#define PE_PARAMETER_CALL8(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8
#define PE_PARAMETER_CALL9(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9
#define PE_PARAMETER_CALL10(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10
#define PE_PARAMETER_CALL11(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11

#define PE_PARAMETER_ADD_CALL0(p)
#define PE_PARAMETER_ADD_CALL1(p) ,p##1
#define PE_PARAMETER_ADD_CALL2(p) ,p##1, p##2
#define PE_PARAMETER_ADD_CALL3(p) ,p##1, p##2, p##3
#define PE_PARAMETER_ADD_CALL4(p) ,p##1, p##2, p##3, p##4
#define PE_PARAMETER_ADD_CALL5(p) ,p##1, p##2, p##3, p##4, p##5
#define PE_PARAMETER_ADD_CALL6(p) ,p##1, p##2, p##3, p##4, p##5, p##6
#define PE_PARAMETER_ADD_CALL7(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7
#define PE_PARAMETER_ADD_CALL8(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8
#define PE_PARAMETER_ADD_CALL9(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9
#define PE_PARAMETER_ADD_CALL10(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10
#define PE_PARAMETER_ADD_CALL11(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11

#define PE_PARAMETER_DEF0(T,p)
#define PE_PARAMETER_DEF1(T,p) T##1 p##1
#define PE_PARAMETER_DEF2(T,p) T##1 p##1, T##2 p##2
#define PE_PARAMETER_DEF3(T,p) T##1 p##1, T##2 p##2, T##3 p##3
#define PE_PARAMETER_DEF4(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4
#define PE_PARAMETER_DEF5(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5
#define PE_PARAMETER_DEF6(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6
#define PE_PARAMETER_DEF7(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7
#define PE_PARAMETER_DEF8(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8
#define PE_PARAMETER_DEF9(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9
#define PE_PARAMETER_DEF10(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10
#define PE_PARAMETER_DEF11(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11

#define PE_MEMBER_DEF0(T,p)
#define PE_MEMBER_DEF1(T,p) T##1 p##1;
#define PE_MEMBER_DEF2(T,p) T##1 p##1; T##2 p##2;
#define PE_MEMBER_DEF3(T,p) T##1 p##1; T##2 p##2; T##3 p##3;
#define PE_MEMBER_DEF4(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4;
#define PE_MEMBER_DEF5(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5;
#define PE_MEMBER_DEF6(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6;
#define PE_MEMBER_DEF7(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7;
#define PE_MEMBER_DEF8(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8;
#define PE_MEMBER_DEF9(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9;
#define PE_MEMBER_DEF10(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9; T##10 p##10;
#define PE_MEMBER_DEF11(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9; T##10 p##10; T##11 p##11;

#define PE_INITIALIZE0(m,p)
#define PE_INITIALIZE1(m,p) : m##1(p##1)
#define PE_INITIALIZE2(m,p) : m##1(p##1), m##2(p##2)
#define PE_INITIALIZE3(m,p) : m##1(p##1), m##2(p##2), m##3(p##3)
#define PE_INITIALIZE4(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4)
#define PE_INITIALIZE5(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5)
#define PE_INITIALIZE6(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6)
#define PE_INITIALIZE7(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7)
#define PE_INITIALIZE8(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8)
#define PE_INITIALIZE9(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9)
#define PE_INITIALIZE10(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9), m##10(p##10)
#define PE_INITIALIZE11(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9), m##10(p##10), m##11(p##11)

#define PE_PARAMETER_CONSTRUCT0(p)
#define PE_PARAMETER_CONSTRUCT1(p) (p##1)
#define PE_PARAMETER_CONSTRUCT2(p) (p##1,p##2)
#define PE_PARAMETER_CONSTRUCT3(p) (p##1,p##2,p##3)
#define PE_PARAMETER_CONSTRUCT4(p) (p##1,p##2,p##3,p##4)
#define PE_PARAMETER_CONSTRUCT5(p) (p##1,p##2,p##3,p##4,p##5)
#define PE_PARAMETER_CONSTRUCT6(p) (p##1,p##2,p##3,p##4,p##5,p##6)
#define PE_PARAMETER_CONSTRUCT7(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7)
#define PE_PARAMETER_CONSTRUCT8(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8)
#define PE_PARAMETER_CONSTRUCT9(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9)
#define PE_PARAMETER_CONSTRUCT10(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p##10)
#define PE_PARAMETER_CONSTRUCT11(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p##10,p##11)



// class to hold the information about what we're executing to pass to the threadmain
template <typename Functor>
class ParallelExecutorClassBase {
 public:
  int numThreads;
  int id;
  CSObject *cs;

  Functor f;
  void *params;
};



// classes to hold a set of parameters - surely some template generalization of this exists somewhere
class ParallelExecutorParameters0 {
 public:
  ParallelExecutorParameters0() { }
};
#define PE_PARAMETER_CLASS(np)						\
    template PE_TEMPLATE_DEF##np					\
	class ParallelExecutorParameters##np {				\
	    public:							\
		ParallelExecutorParameters##np(PE_PARAMETER_DEF##np(T,&_p))	\
		    PE_INITIALIZE##np(p,_p) { }				\
		PE_MEMBER_DEF##np(T,&p)					\
	};


// functions of type void* f(void*) - these are the threadmain's
#define PE_STUB(np)							\
    template <typename Functor PE_TEMPLATE_ADD_DEF##np>			\
	void* ParallelExecutorStub##np(void *arg) {			\
	ParallelExecutorClassBase<Functor> *pec = (ParallelExecutorClassBase<Functor>*)arg; \
	ParallelExecutorParameters##np PE_TEMPLATE_CALL##np *params = (ParallelExecutorParameters##np PE_TEMPLATE_CALL##np*)pec->params; \
	pec->f(pec->numThreads, pec->id, *(pec->cs) PE_PARAMETER_ADD_CALL##np(params->p)); \
	return 0;							\
    }




// make sure you don't try to recursively call in parallel - I've already made this mistake once!
inline void ParallelExecutorCheck(bool start) {
  static int count = 0;
  static CSObject cs;

  cs.enter();
  if (start) {
    if (count > 0)
      cerr<<"starting parallel executor recursively?!?!"<<endl;
    count++;
  } else {
    count--;
  }
  cs.leave();
}


// the real work - create the set of threads, have them run, then join them back together
#define PARALLEL_EXECUTOR(np)											\
    template <typename Functor PE_TEMPLATE_ADD_DEF##np>					\
	void ParallelExecutor(const int reqThreads, Functor f PE_PARAMETER_ADD_DEF##np(T,&p)) { \
		ParallelExecutorCheck(true);									\
	ParallelExecutorParameters##np PE_TEMPLATE_CALL##np params PE_PARAMETER_CONSTRUCT##np(p); \
	vector< ParallelExecutorClassBase<Functor> > pec(reqThreads);		\
	CSObject cs;														\
	for (int i=0; i<reqThreads; i++) {									\
	    pec[i].numThreads = reqThreads;									\
	    pec[i].id = i;													\
		pec[i].cs = &cs;												\
	    pec[i].f = f;													\
	    pec[i].params = &params;										\
	}																	\
	vector<thlib::Thread *> threads(reqThreads-1);						\
	for (int i=0; i<reqThreads-1; i++) {								\
	    threads[i] = new thlib::Thread(&ParallelExecutorStub##np<Functor PE_TEMPLATE_ADD_CALL##np>, (void*)&pec[i], 0); \
	}																	\
	ParallelExecutorStub##np<Functor PE_TEMPLATE_ADD_CALL##np>(&pec[reqThreads-1]); \
	for (int i=0; i<reqThreads-1; i++) {								\
	    int *status;													\
	    threads[i]->join((void**)&status);								\
	    delete threads[i];												\
	}																	\
	ParallelExecutorCheck(false);										\
    }


// instantiate them for up to 5 parameters

PE_PARAMETER_CLASS(1)
PE_PARAMETER_CLASS(2)
PE_PARAMETER_CLASS(3)
PE_PARAMETER_CLASS(4)
PE_PARAMETER_CLASS(5)
PE_PARAMETER_CLASS(6)
PE_PARAMETER_CLASS(7)
PE_PARAMETER_CLASS(8)
PE_PARAMETER_CLASS(9)
PE_PARAMETER_CLASS(10)
PE_PARAMETER_CLASS(11)

PE_STUB(0)
PE_STUB(1)
PE_STUB(2)
PE_STUB(3)
PE_STUB(4)
PE_STUB(5)
PE_STUB(6)
PE_STUB(7)
PE_STUB(8)
PE_STUB(9)
PE_STUB(10)
PE_STUB(11)

PARALLEL_EXECUTOR(0)
PARALLEL_EXECUTOR(1)
PARALLEL_EXECUTOR(2)
PARALLEL_EXECUTOR(3)
PARALLEL_EXECUTOR(4)
PARALLEL_EXECUTOR(5)
PARALLEL_EXECUTOR(6)
PARALLEL_EXECUTOR(7)
PARALLEL_EXECUTOR(8)
PARALLEL_EXECUTOR(9)
PARALLEL_EXECUTOR(10)
PARALLEL_EXECUTOR(11)






// I think the std::mem_fun stuff is similar, but doesn't take any number of parameters.
// to use with the parallel executor functions, you need instantiation for 2 extra
// variables since the thread id and num threads parameters are not treated as special.
template <typename Class, typename Functor>
class ClassFunctor {
 public:
  ClassFunctor() { }
  ClassFunctor(Class *_c, Functor _f) : c(_c), f(_f) { }
  
  void operator()() { (c->*f)(); }


#define PE_CLASS_FUNCTOR(np)						\
    template PE_TEMPLATE_DEF##np					\
	void operator() (PE_PARAMETER_DEF##np(T,&p)) {			\
	(c->*f)(PE_PARAMETER_CALL##np(p));				\
    }									\

PE_CLASS_FUNCTOR(1)
PE_CLASS_FUNCTOR(2)
PE_CLASS_FUNCTOR(3)
PE_CLASS_FUNCTOR(4)
PE_CLASS_FUNCTOR(5)
PE_CLASS_FUNCTOR(6)
PE_CLASS_FUNCTOR(7)
PE_CLASS_FUNCTOR(8)
PE_CLASS_FUNCTOR(9)
PE_CLASS_FUNCTOR(10)
PE_CLASS_FUNCTOR(11)


 private:
  Class *c;
  Functor f;
};


// this will autodetect/deduce the types, so it's easy to place inline (create a temp object) without
// having to figure out what the damn syntax is for member function pointers
template <typename Class, typename Functor>
ClassFunctor<Class, Functor> makeClassFunctor(Class *c, Functor f) {
  return ClassFunctor<Class,Functor>(c,f);
}


extern int idealNumThreads;



#endif

