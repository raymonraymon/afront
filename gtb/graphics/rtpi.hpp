
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


#ifndef GTB_RTPI_INCLUDED
#define GTB_RTPI_INCLUDED

#include <gtb/common.hpp>
#include <gtb/real/real.hpp>
#include <gtb/graphics/point3.hpp>


GTB_BEGIN_NAMESPACE


// Range, theta, phi, intensity
class Rtpi {
public:
	Rtpi();
	Rtpi(real_t r, real_t t, real_t p, int i);

	real_t r() const;
	real_t t() const;
	real_t p() const;
	int i() const;

	void set_r(real_t);
	void set_t(real_t);
	void set_p(real_t);
	void set_i(int);
	void reset(real_t r, real_t t, real_t p, int i);

	Point3 point() const;

	void read(FILE *fp);
	void write(FILE *fp) const;

protected:
	real_t _r, _t, _p;
	int _i;
};


GTB_END_NAMESPACE

#ifndef OUTLINE
#include <gtb/graphics/rtpi.ipp>
#endif

#endif // GTB_RTPI_INCLUDED
