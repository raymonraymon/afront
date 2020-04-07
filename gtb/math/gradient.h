
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



#ifndef GRADIENT_H
#define GRADIENT_H


template <typename T, int N>
class Gradient {
public:

	T value;
	T gradient[N];


	Gradient() { }
	Gradient(const T &v) { 
		value = v;
		for (int i=0; i<N; i++) { gradient[i] = 0.0; }
	}
      

	inline Gradient<T,N> operator+(const Gradient<T,N> &rhs) const {
		Gradient<T,N> ret;
		ret.value = value + rhs.value;
		for (int i=0; i<N; i++) { ret.gradient[i] = gradient[i] + rhs.gradient[i]; }
		return ret;
	}
	inline Gradient<T,N>& operator+=(const Gradient<T,N> &rhs) {
		value += rhs.value;
		for (int i=0; i<N; i++) { gradient[i] += rhs.gradient[i]; }
		return *this;
	}

	inline Gradient<T,N> operator-(const Gradient<T,N> &rhs) const {
		Gradient<T,N> ret;
		ret.value = value - rhs.value;
		for (int i=0; i<N; i++) { ret.gradient[i] = gradient[i] - rhs.gradient[i]; }
		return ret;
	}
	inline Gradient<T,N>& operator-=(const Gradient<T,N> &rhs) {
		value -= rhs.value;
		for (int i=0; i<N; i++) { gradient[i] -= rhs.gradient[i]; }
		return *this;
	}

	inline Gradient<T,N> operator*(const Gradient<T,N> &rhs) const {
		// product rule
		Gradient<T,N> ret;
		for (int i=0; i<N; i++) { ret.gradient[i] = rhs.value*gradient[i] + value*rhs.gradient[i]; }
		ret.value = value * rhs.value;
		return ret;
	}
	inline Gradient<T,N>& operator*=(const Gradient<T,N> &rhs) {
		// product rule
		for (int i=0; i<N; i++) { gradient[i] = rhs.value*gradient[i] + value*rhs.gradient[i]; }
		value = value * rhs.value;
		return *this;
	}

	inline Gradient<T,N> operator/(const Gradient<T,N> &rhs) const {
		// quotient rule
		Gradient<T,N> ret;
		for (int i=0; i<N; i++) { ret.gradient[i] = (rhs.value*gradient[i] - value*rhs.gradient[i]) / (rhs.value*rhs.value); }
		ret.value = value / rhs.value;
		return ret;
	}
	inline Gradient<T,N>& operator/=(const Gradient<T,N> &rhs) {
		// quotient rule
		for (int i=0; i<N; i++) { gradient[i] = (rhs.value*gradient[i] - value*rhs.gradient[i]) / (rhs.value*rhs.value); }
		value = value / rhs.value;
		return *this;
	}


	inline Gradient<T,N> ChainRule(const T &fx, const T &fpx) const {
		Gradient<T,N> ret;
		ret.value = fx;
		for (int i=0; i<N; i++) { ret.gradient[i] = fpx*gradient[i]; }
		return ret;
	}


	// negate
	inline Gradient<T,N> operator-() const {
		Gradient<T,N> ret;
		ret.value = -value;
		for (int i=0; i<N; i++) { ret.gradient[i] = -gradient[i]; }
		return ret;
	}




	//
	// ops on constant builtin types
	//


#define GRADIENT_OPERATOR_EQ(type) \
	inline Gradient<T,N>& operator=(type c) {\
		value = c; \
		for (int i=0; i<N; i++) { gradient[i] = 0.0; } \
		return *this; \
	}

	// multiplication by constants
#define GRADIENT_OPERATOR_MULT(type) \
	inline Gradient<T,N> operator*(type c) const { \
		Gradient<T,N> ret; \
		for (int i=0; i<N; i++) { ret.gradient[i] = gradient[i] * c; } \
		ret.value = value * c; \
		return ret; \
	}

#define GRADIENT_OPERATOR_MULTEQ(type) \
	inline Gradient<T,N> operator*=(type c) { \
		for (int i=0; i<N; i++) { gradient[i] = gradient[i] * c; } \
		value = value * c; \
		return *this; \
	}

#define GRADIENT_FRIEND_MULT(type) \
	inline friend Gradient<T,N> operator*(type c, const Gradient<T,N> &g) { \
		return (g*c); \
	}


	// division by constants
#define GRADIENT_OPERATOR_DIV(type) \
	inline Gradient<T,N> operator/(type c) const { \
		Gradient<T,N> ret; \
		for (int i=0; i<N; i++) { ret.gradient[i] = gradient[i] / c; } \
		ret.value = value / c; \
		return ret; \
	}

#define GRADIENT_OPERATOR_DIVEQ(type) \
	inline Gradient<T,N> operator/=(type c) { \
		for (int i=0; i<N; i++) { gradient[i] = gradient[i] / c; } \
		value = value / c; \
		return *this; \
	}

#define GRADIENT_FRIEND_DIV(type) \
	inline friend Gradient<T,N> operator/(type c, const Gradient<T,N> &g) { \
		return (c*inverse(g)); \
	}


	// addition of constants
#define GRADIENT_OPERATOR_PLUS(type) \
	inline Gradient<T,N> operator+(type c) const { \
		Gradient<T,N> ret = *this; \
		ret.value = value + c; \
		return ret; \
	}

#define GRADIENT_OPERATOR_PLUSEQ(type) \
	inline Gradient<T,N> operator+=(type c) const { \
		value = value + c; \
		return *this; \
	}

#define GRADIENT_FRIEND_PLUS(type) \
	inline friend Gradient<T,N> operator+(type c, const Gradient<T,N> &g) { \
		return (g+c); \
	}



	// subtraction of constants
#define GRADIENT_OPERATOR_SUB(type) \
	inline Gradient<T,N> operator-(type c) const { \
		Gradient<T,N> ret = *this; \
		ret.value = value - c; \
		return ret; \
	}

#define GRADIENT_OPERATOR_SUBEQ(type) \
	inline Gradient<T,N> operator-=(type c) const { \
		value = value - c; \
		return *this; \
	}

#define GRADIENT_FRIEND_SUB(type) \
	inline friend Gradient<T,N> operator-(type c, const Gradient<T,N> &g) { \
		Gradient<T,N> ret; \
		ret.value = c - g.value; \
		for (int i=0; i<N; i++) { ret.gradient[i] = -g.gradient[i]; } \
		return ret; \
	}


#define GRADIENT_DEFINE_CONST_OPS(type) \
	GRADIENT_OPERATOR_EQ(type) \
	GRADIENT_OPERATOR_MULT(type) \
	GRADIENT_OPERATOR_MULTEQ(type) \
	GRADIENT_FRIEND_MULT(type) \
	GRADIENT_OPERATOR_DIV(type) \
	GRADIENT_OPERATOR_DIVEQ(type) \
	GRADIENT_FRIEND_DIV(type) \
	GRADIENT_OPERATOR_PLUS(type) \
	GRADIENT_OPERATOR_PLUSEQ(type) \
	GRADIENT_FRIEND_PLUS(type) \
	GRADIENT_OPERATOR_SUB(type) \
	GRADIENT_OPERATOR_SUBEQ(type) \
	GRADIENT_FRIEND_SUB(type)


	GRADIENT_DEFINE_CONST_OPS(int)
	GRADIENT_DEFINE_CONST_OPS(float)
	GRADIENT_DEFINE_CONST_OPS(double)

#ifndef WIN32
	GRADIENT_DEFINE_CONST_OPS(long double)
#endif


};


template <typename T, int N>	// much faster and stabler than 1/x or pow(x,-1) ?
inline Gradient<T,N> inverse(const Gradient<T,N> &v) { T inverse=1.0/v.value;  return v.ChainRule(inverse, -inverse*inverse); }

template <typename T, int N>	// much faster and stabler than pow(x,1/2) ?
inline Gradient<T,N> sqrt(const Gradient<T,N> &v) { T s=sqrt(v.value); return v.ChainRule(s, 0.5/s); }



template <typename T, int N>
inline Gradient<T,N> sin(const Gradient<T,N> &v) { return v.ChainRule(sin(v.value),  cos(v.value)); }

template <typename T, int N>
inline Gradient<T,N> cos(const Gradient<T,N> &v) { return v.ChainRule(cos(v.value), -sin(v.value)); }

template <typename T, int N>
inline Gradient<T,N> tan(const Gradient<T,N> &v) { T tangent=tan(v.value); return v.ChainRule(tangent, 1+tangent*tangent); }


// hack hack hack
inline double pow(double b, float e) { return pow((double)b, (double)e); }
inline double pow(float b, double e) { return pow((double)b, (double)e); }
template <typename T, int N>
inline Gradient<T,N> pow(const Gradient<T,N> &v, double e) { return v.ChainRule(pow(v.value, e), e*pow(v.value, e-1)); }

template <typename T, int N>
inline Gradient<T,N> exp(const Gradient<T,N> &v) { T e=exp(v.value); return v.ChainRule(e, e); }



// comparisons
template <typename T, int N>
inline bool operator<(const Gradient<T,N> &l, const Gradient<T,N> &r) { return (l.value<r.value); }
template <typename T, int N>
inline bool operator<(float l, const Gradient<T,N> &r) { return (l<r.value); }
template <typename T, int N>
inline bool operator<(const Gradient<T,N> &l, float r) { return (l.value<r); }

template <typename T, int N>
inline bool operator<=(const Gradient<T,N> &l, const Gradient<T,N> &r) { return (l.value<=r.value); }
template <typename T, int N>
inline bool operator<=(float l, const Gradient<T,N> &r) { return (l<=r.value); }
template <typename T, int N>
inline bool operator<=(const Gradient<T,N> &l, float r) { return (l.value<=r); }

template <typename T, int N>
inline bool operator>(const Gradient<T,N> &l, const Gradient<T,N> &r) { return (l.value>r.value); }
template <typename T, int N>
inline bool operator>(float l, const Gradient<T,N> &r) { return (l>r.value); }
template <typename T, int N>
inline bool operator>(const Gradient<T,N> &l, float r) { return (l.value>r); }

template <typename T, int N>
inline bool operator>=(const Gradient<T,N> &l, const Gradient<T,N> &r) { return (l.value>=r.value); }
template <typename T, int N>
inline bool operator>=(float l, const Gradient<T,N> &r) { return (l>=r.value); }
template <typename T, int N>
inline bool operator>=(const Gradient<T,N> &l, float r) { return (l.value>=r); }

template <typename T, int N>
inline bool operator==(const Gradient<T,N> &l, const Gradient<T,N> &r) { return (l.value==r.value); }
template <typename T, int N>
inline bool operator==(float l, const Gradient<T,N> &r) { return (l==r.value); }
template <typename T, int N>
inline bool operator==(const Gradient<T,N> &l, float r) { return (l.value==r); }



template<typename T>
T fabs(T x) { return (x<0) ? -x : x; }


#endif





















