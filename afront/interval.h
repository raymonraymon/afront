
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



#ifndef INTERVAL_H
#define INTERVAL_H


// very simple interval arithmetic class - doesn't change any rounding modes at all!
class Interval {
public:
	real_type min, max;

	Interval() {}
	Interval(const real_type &x) {
		min = max = x;
	}

	Interval(const real_type &min, const real_type &max) {
		this->min = min;
		this->max = max;
	}
	
	inline real_type& operator[](int i) {
		if (i==0)
			return min;
		if (i==1)
			return max;
		cerr<<"bad interval index!!"<<endl;
		exit(0);
	}

	inline const real_type& operator[](int i) const {
		if (i==0)
			return min;
		if (i==1)
			return max;
		cerr<<"bad interval index!!"<<endl;
		exit(0);
	}

	inline Interval& operator=(int c) {
		min = max = (real_type)c;
		return *this;
	}
	inline Interval& operator=(float c) {
		min = max = (real_type)c;
		return *this;
	}
	inline Interval& operator=(double c) {
		min = max = (real_type)c;
		return *this;
	}
	inline Interval& operator=(const Interval &i) {
		min = i.min;
		max = i.max;
		return *this;
	}


	// addition
	template <typename T>
	inline Interval operator+(const T &c) const {
		return Interval(min+c, max+c);
	}
	inline Interval operator+(const Interval &rhs) const {
		return Interval(min+rhs.min, max+rhs.max);
	}

	template <typename T>
	inline Interval& operator+=(const T &c) {
		min += c;
		max += c;
		return *this;
	}
	inline Interval& operator+=(const Interval &rhs) {
		min += rhs.min;
		max += rhs.max;
		return *this;
	}

	template <typename T>
	inline friend Interval operator+(const T &c, const Interval &i) {
		return (i+c);
	}


	// subtraction
	template <typename T>
	inline Interval operator-(const T &c) const {
		return Interval(min-c, max-c);
	}
	inline Interval operator-(const Interval &rhs) const {
		return Interval(min-rhs.min, max-rhs.max);
	}

	template <typename T>
	inline Interval& operator-=(const T &c) {
		min -= c;
		max -= c;
		return *this;
	}
	inline Interval& operator-=(const Interval &rhs) {
		min -= rhs.min;
		max -= rhs.max;
		return *this;
	}

	template <typename T>
	inline friend Interval operator-(const T &c, const Interval &i) {
		return (i-c);
	}


	// multiplication
	template <typename T>
	inline Interval operator*(const T &c) const {
		if (c<0)
			return Interval(max*c, min*c);
		else
			return Interval(min*c, max*c);
	}
	inline Interval operator*(const Interval &rhs) const {
		return Interval(std::min(std::min(min*rhs.min, min*rhs.max),
								 std::min(max*rhs.min, max*rhs.max)),
						std::max(std::max(min*rhs.min, min*rhs.max),
								 std::max(max*rhs.min, max*rhs.max)));

	}
	template <typename T>
	inline Interval& operator*=(const T &c) {
		if (c<0)
			std::swap(min,max);
		min *= c;
		max *= c;
		return *this;
	}
	inline Interval& operator*=(const Interval &rhs) {
		real_type nmin = std::min(std::min(min*rhs.min, min*rhs.max),
								  std::min(max*rhs.min, max*rhs.max));
		real_type nmax = std::max(std::max(min*rhs.min, min*rhs.max),
								  std::max(max*rhs.min, max*rhs.max));

		min = nmin;
		max = nmax;
		return *this;
	}

	template <typename T>
	inline friend Interval operator*(const T &c, const Interval &i) {
		return (i*c);
	}




	// division
	template <typename T>
	inline Interval operator/(const T &c) const {
		if (c==0) {
			return Interval(NAN, NAN);
		}
		else if (c<0)
			return Interval(max/c, min/c);
		else
			return Interval(min/c, max/c);
	}
	inline Interval operator/(const Interval &rhs) const {
		if (rhs.min<0 || rhs.max>0) {
			// no divide by 0!
			return Interval(std::min(std::min(min/rhs.min, min/rhs.max),
									 std::min(max/rhs.min, max/rhs.max)),
							std::max(std::max(min/rhs.min, min/rhs.max),
									 std::max(max/rhs.min, max/rhs.max)));
		} else {
			// trying to divide by 0
			return Interval(NAN, NAN);
		}
	}
	template <typename T>
	inline Interval& operator/=(const T &c) {
		if (c==0) {
			min = max = NAN;
		}
		else if (c<0) {
			std::swap(min,max);
			min /= c;
			max /= c;
		} else {
			min /= c;
			max /= c;
		}
		return *this;
	}
	inline Interval& operator/=(const Interval &rhs) {
		if (rhs.min<0 || rhs.max>0) {
			// no divide by 0
			real_type nmin = std::min(std::min(min/rhs.min, min/rhs.max),
									  std::min(max/rhs.min, max/rhs.max));
			real_type nmax = std::max(std::max(min/rhs.min, min/rhs.max),
									  std::max(max/rhs.min, max/rhs.max));
			
			min = nmin;
			max = nmax;
		} else {
			// trying to divide by 0
			min = NAN;
			max = NAN;
		}
		return *this;
	}
	template <typename T>
	inline friend Interval operator/(const T &c, const Interval &i) {
		return (i/c);
	}



	static Interval Intersect(const Interval &x, const Interval &y) {
		return Interval(std::max(x.min, y.min), std::min(x.max, y.max));
	}

	static Interval Union(const Interval &x, const Interval &y) {
		return Interval(std::min(x.min, y.min), std::max(x.max, y.max));
	}
};


// specialize square and cube
inline Interval square(const Interval &x) {
	if (x.max < 0) {
		// all negative
		return Interval(square(x.max), square(x.min));
	} else if (x.min>0) {
		// all positive
		return Interval(square(x.min), square(x.max));
	} else {
		// spans 0
		return Interval(0, std::max(square(x.min), square(x.max)));
	}
}

inline Interval cube(const Interval &x) {
	// easy - monotonic increasing
	return Interval(cube(x.min), cube(x.max));
}


inline Interval quadratic(const real_type c[3], const Interval &x) {

	Interval ret(quadratic(c, x.min), quadratic(c,x.max));
	if (ret.min > ret.max)
		std::swap(ret.min, ret.max);

	// check if we have a min or max to deal with
	if (c[0] != 0) {
		real_type xmin = -c[1] / (2*c[0]);

		if (x.min < xmin && x.max>xmin) {
			real_type f = quadratic(c, xmin);
			if (f<ret.min)
				ret.min = f;
			else if (f>ret.max)
				ret.max = f;
		}
	}

	return ret;
}


template <typename T1>
inline Interval cubic(const T1 c[4], const Interval &x) {

	if (c[0] == 0)
		return quadratic(&c[1], x);

	Interval ret(cubic(c, x.min), cubic(c,x.max));
	if (ret.min > ret.max)
		std::swap(ret.min, ret.max);


	real_type cderiv[3] = { 3*c[0], 2*c[1], c[2] };

	real_type descrim = square(cderiv[1]) - 4*cderiv[0]*cderiv[2];

	if (descrim > 0) {

		// we have actual root in the derivative - consider the local min/max of the cubic
		real_type t;
		if (cderiv[1]>0) {
			t = -(cderiv[1] + sqrt(descrim)) / 2;
		} else {
			t = -(cderiv[1] - sqrt(descrim)) / 2;
		}

		real_type exts[2] = { t/cderiv[0], cderiv[2]/t };

		for (int i=0; i<2; i++) {
			if (x.min < exts[i] && x.max > exts[i]) {
				real_type f = cubic(c, exts[i]);
				if (f<ret.min)
					ret.min = f;
				else if (f>ret.max)
					ret.max = f;
			}
		}
	}

	return ret;
}




// sqrt is monotonically increasing so it's easy
inline Interval sqrt(const Interval &i) {
	if (i.min<0) {
		return Interval(NAN, NAN);
	}

	return Interval(sqrt(i.min), sqrt(i.max));
}


inline Interval abs(const Interval &i) {
	if (i.min > 0)
		return Interval(i.min, i.max);

	if (i.max < 0)
		return Interval(-i.max, -i.min);

	// cross 0
	return Interval(0, std::max(-i.min, i.max));
}



inline Interval pow(const Interval &x, int p) {

	if (p==0)
		return Interval(1,1);

	// if its odd it's monotonically increasing
	if (p&1) {
		return Interval(pow(x.min, p), pow(x.max, p));
	} else {
		if (x.min>0)
			return Interval(pow(x.min, p), pow(x.max, p)); // monotonically increasing

		if (x.max<0)
			return Interval(pow(x.max, p), pow(x.min, p)); // monotonically decreasing

		// otherwise, it has a min at 0
		return Interval(0, std::max(pow(x.min, p), pow(x.max, p)));
	}
}



inline int interval_sign(const Interval &i) {
	if (i.max<0)
		return -1;
	if (i.min>0)
		return 1;
	return 0;
}





#endif
