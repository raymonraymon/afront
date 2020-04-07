
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
#ifndef WIN32
#include <gtb/math/math.hpp>
#include <gtb/math/matrix4.hpp>
#include <iostream>
#endif // WIN32

using namespace gtb;
using namespace std;


int main(int, char *[])
{
	cout << "identity =\n" << MATRIX4_IDENTITY << "\n";

	Matrix4 m1;
	cout << "m1 =\n" << m1 << "\n";

	real_t a2[4][4] = {
		{200, 201, 202, 203},
		{210, 211, 212, 213},
		{220, 221, 222, 223},
		{230, 231, 232, 233}
	};
	Matrix4 m2(a2);
	cout << "m2 =\n" << m2 << "\n";

	real_t a3[16] = {
		300, 301, 302, 303,
		310, 311, 312, 313,
		320, 321, 322, 323,
		330, 331, 332, 333
	};
	Matrix4 m3(a3);
	cout << "m3 =\n" << m3 << "\n";

	Matrix4 m4(400, 401, 402, 403,
		   410, 411, 412, 413,
		   420, 421, 422, 423,
		   430, 431, 432, 433);
	cout << "m4 =\n" << m4 << "\n";

	Matrix4 m5;
	cout << "m5 =\n" << m5 << "\n";
	m5 = m4;
	cout << "m5 =\n" << m5 << "\n";

	Matrix4 m6(m4);
	m6.make_identity();
	cout << "m6 =\n" << m6 << "\n";

	Matrix4 m7(m4);
	m7.negate();
	cout << "m7 =\n" << m7 << "\n";

	Matrix4 m8 = m4.transpose();
	cout << "m8 =\n" << m8 << "\n";

	real_t a9[16] = { 
		1, 6, 2, 4,
		3, 19, 4, 15,
		1, 4, 8, 12,
		5, 33, 9, 3
	};
	Matrix4 m9(a9);
	cout << "m9 =\n" << m9 << "\n";

	Matrix4 m10 = m9.inverse();
	cout << "m10 =\n" << m10 << "\n";
	cout << "m9 * m10 =\n" << m9 * m10 << "\n";
	cout << "det(m9) = " << m9.det() << "\n";

	real::set_eps(1.0e-5);
	cout << "(m9 * m10 == identity) = ";
	cout << (m9 * m10 == MATRIX4_IDENTITY) << "\n";

	Matrix4 m11;
	cout << "(m11 == zero) = ";
	cout << (m11 == MATRIX4_ZERO) << "\n";

	cout << "-m4 = \n" << -m4 << "\n";
	cout << "m4 += m4 \n" << (m4 += m4) << "\n";
	cout << "m4 -= 0.5 * m4 \n" << (m4 -= 0.5 * m4) << "\n";
	cout << "m4 *= 2 =\n" << (m4 *= 2) << "\n";
	cout << "m4 *= identity =\n" << (m4 *= MATRIX4_IDENTITY) << "\n";
	cout << "m4 *= 0.5 =\n" << (m4 *= 0.5) << "\n";

	cout << "m4 + m4 =\n" << m4 + m4 << "\n";
	cout << "m4 - m4 =\n" << m4 - m4 << "\n";
	cout << "m4 * m4 =\n" << m4 * m4 << "\n";
	cout << "2 * m4 =\n" << 2 * m4 << "\n";
	cout << "m4 * 2 =\n" << m4 * 2 << "\n";

	return EXIT_SUCCESS;
}
