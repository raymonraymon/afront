
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
#include <gtb/math/matrix4.hpp>
#include <gtb/error/error.hpp>
#include <gtb/math/math.hpp>
#include <math.h>
//#include <stdlib.h>
#endif // WIN32

#ifdef OUTLINE
#define inline
#include <gtb/math/matrix4.ipp>
#undef inline
#endif

using namespace std;

GTB_BEGIN_NAMESPACE


template<class T>
tMatrix4<T>::tMatrix4()
{
	_m[0][0] = 0.0;
	_m[0][1] = 0.0;
	_m[0][2] = 0.0;
	_m[0][3] = 0.0;

	_m[1][0] = 0.0;
	_m[1][1] = 0.0;
	_m[1][2] = 0.0;
	_m[1][3] = 0.0;

	_m[2][0] = 0.0;
	_m[2][1] = 0.0;
	_m[2][2] = 0.0;
	_m[2][3] = 0.0;

	_m[3][0] = 0.0;
	_m[3][1] = 0.0;
	_m[3][2] = 0.0;
	_m[3][3] = 0.0;
}


template<class T>
tMatrix4<T>::tMatrix4(value_type m[4][4])
{
	_m[0][0] = m[0][0];
	_m[0][1] = m[0][1];
	_m[0][2] = m[0][2];
	_m[0][3] = m[0][3];

	_m[1][0] = m[1][0];
	_m[1][1] = m[1][1];
	_m[1][2] = m[1][2];
	_m[1][3] = m[1][3];

	_m[2][0] = m[2][0];
	_m[2][1] = m[2][1];
	_m[2][2] = m[2][2];
	_m[2][3] = m[2][3];

	_m[3][0] = m[3][0];
	_m[3][1] = m[3][1];
	_m[3][2] = m[3][2];
	_m[3][3] = m[3][3];
}

template<class T>
tMatrix4<T>::tMatrix4(value_type m[16])
{
	_m[0][0] = m[0];
	_m[0][1] = m[1];
	_m[0][2] = m[2];
	_m[0][3] = m[3];

	_m[1][0] = m[4];
	_m[1][1] = m[5];
	_m[1][2] = m[6];
	_m[1][3] = m[7];

	_m[2][0] = m[8];
	_m[2][1] = m[9];
	_m[2][2] = m[10];
	_m[2][3] = m[11];

	_m[3][0] = m[12];
	_m[3][1] = m[13];
	_m[3][2] = m[14];
	_m[3][3] = m[15];
}

template<class T>
tMatrix4<T>::tMatrix4(value_type m00, value_type m01, value_type m02, value_type m03,
		 value_type m10, value_type m11, value_type m12, value_type m13,
		 value_type m20, value_type m21, value_type m22, value_type m23,
		 value_type m30, value_type m31, value_type m32, value_type m33)
{
	_m[0][0] = m00;
	_m[0][1] = m01;
	_m[0][2] = m02;
	_m[0][3] = m03;

	_m[1][0] = m10;
	_m[1][1] = m11;
	_m[1][2] = m12;
	_m[1][3] = m13;

	_m[2][0] = m20;
	_m[2][1] = m21;
	_m[2][2] = m22;
	_m[2][3] = m23;

	_m[3][0] = m30;
	_m[3][1] = m31;
	_m[3][2] = m32;
	_m[3][3] = m33;
}

template<class T>
tMatrix4<T>::tMatrix4(const tMatrix4 &m)
{
	_m[0][0] = m._m[0][0];
	_m[0][1] = m._m[0][1];
	_m[0][2] = m._m[0][2];
	_m[0][3] = m._m[0][3];

	_m[1][0] = m._m[1][0];
	_m[1][1] = m._m[1][1];
	_m[1][2] = m._m[1][2];
	_m[1][3] = m._m[1][3];

	_m[2][0] = m._m[2][0];
	_m[2][1] = m._m[2][1];
	_m[2][2] = m._m[2][2];
	_m[2][3] = m._m[2][3];

	_m[3][0] = m._m[3][0];
	_m[3][1] = m._m[3][1];
	_m[3][2] = m._m[3][2];
	_m[3][3] = m._m[3][3];
}

template<class T>
tMatrix4<T> &tMatrix4<T>::operator=(const tMatrix4 &m)
{
	if (&m != this) {
		_m[0][0] = m._m[0][0];
		_m[0][1] = m._m[0][1];
		_m[0][2] = m._m[0][2];
		_m[0][3] = m._m[0][3];

		_m[1][0] = m._m[1][0];
		_m[1][1] = m._m[1][1];
		_m[1][2] = m._m[1][2];
		_m[1][3] = m._m[1][3];

		_m[2][0] = m._m[2][0];
		_m[2][1] = m._m[2][1];
		_m[2][2] = m._m[2][2];
		_m[2][3] = m._m[2][3];

		_m[3][0] = m._m[3][0];
		_m[3][1] = m._m[3][1];
		_m[3][2] = m._m[3][2];
		_m[3][3] = m._m[3][3];
	}
	return *this;
}

template<class T>
tMatrix4<T> &tMatrix4<T>::make_identity()
{
	_m[0][0] = 1.0;
	_m[0][1] = 0.0;
	_m[0][2] = 0.0;
	_m[0][3] = 0.0;

	_m[1][0] = 0.0;
	_m[1][1] = 1.0;
	_m[1][2] = 0.0;
	_m[1][3] = 0.0;

	_m[2][0] = 0.0;
	_m[2][1] = 0.0;
	_m[2][2] = 1.0;
	_m[2][3] = 0.0;

	_m[3][0] = 0.0;
	_m[3][1] = 0.0;
	_m[3][2] = 0.0;
	_m[3][3] = 1.0;

	return *this;
}

template<class T>
tMatrix4<T> &tMatrix4<T>::negate()
{
	_m[0][0] = -_m[0][0];
	_m[0][1] = -_m[0][1];
	_m[0][2] = -_m[0][2];
	_m[0][3] = -_m[0][3];

	_m[1][0] = -_m[1][0];
	_m[1][1] = -_m[1][1];
	_m[1][2] = -_m[1][2];
	_m[1][3] = -_m[1][3];

	_m[2][0] = -_m[2][0];
	_m[2][1] = -_m[2][1];
	_m[2][2] = -_m[2][2];
	_m[2][3] = -_m[2][3];

	_m[3][0] = -_m[3][0];
	_m[3][1] = -_m[3][1];
	_m[3][2] = -_m[3][2];
	_m[3][3] = -_m[3][3];

	return *this;
}

template<class T>
bool tMatrix4<T>::LU_decomposition(int p[4])
{
	int i, j, k;

	for (j = 0; j < 3; j++) {

		// Find line of pivot
		p[j] = j;
		for (i = j + 1; i < 4; i++) {
			if (fabs(_m[i][j]) > fabs(_m[p[j]][j])) {
				p[j] = i;
			}
		}

		// Swap lines if necessary
		if (p[j] != j) {
			i = p[j];
			for (k = 0; k < 4; k++) {
				value_type t = _m[i][k];
				_m[i][k] = _m[j][k];
				_m[j][k] = t;
			}
		}

		// Check if matrix is singular
		if (real::is_zero(_m[j][j])) {
			GTB_WARNING("singular matrix");
			return false;
		}

		// Eliminate elements below diagonal
		for (i = j + 1; i < 4; i++) {
			_m[i][j] = -_m[i][j] / _m[j][j];
			for (k = j + 1; k < 4; k++) {
				_m[i][k] += _m[i][j] * _m[j][k];
			}
		}
	}
	return true;
}

template<class T>
void tMatrix4<T>::LU_back_substitution(const int p[4],
				   tMatrix4 &b,
				   tMatrix4 &x)
{
	int i, j, k;

	// Swap lines of b when necessary
	for (i = 0; i < 3; i++) {
		if (p[i] != i) {
			for (j = 0, k = p[i]; j < 4; j++) {
				value_type t = b._m[i][j];
				b._m[i][j] = b._m[k][j];
				b._m[k][j] = t;
			}
		}
	}

	// Apply multipliers to b
	for (j = 0; j < 3 ; j++) {
		for (i = j + 1; i < 4; i++) {
			for (k = 0; k < 4; k++) {
				b._m[i][k] += _m[i][j] * b._m[j][k];
			}
		}
	}

	// Back substitution
	for (k = 0; k < 4; k++) {
		for (i = 3; i >= 0; i--) {
			x._m[i][k] = b._m[i][k];
			for (j = i + 1; j < 4; j++) {
				x._m[i][k] -= _m[i][j] * x._m[j][k];
			}
			assert(!real::is_zero(_m[i][i]));
			x._m[i][k] /= _m[i][i];
		}
	}
}

template<class T>
tMatrix4<T> tMatrix4<T>::inverse() const
{
	tMatrix4 a(*this);
	tMatrix4 b(MATRIX4_IDENTITY);
	tMatrix4 x;
	int p[4];

	a.LU_decomposition(p);
	a.LU_back_substitution(p, b, x);
	return x;
}

template<class T>
typename tMatrix4<T>::value_type tMatrix4<T>::det() const
{
	value_type d;
	tMatrix4 t(*this);
	int p[4];

	t.LU_decomposition(p);
	d = 1.0;
	for (int i = 0; i < 3; i++) {
		d *= t._m[i][i];
		if (p[i] != i) {
			d = -d;
		}
	}
	d *= t._m[3][3];
	return d;
}

template<class T>
void tMatrix4<T>::print(FILE *fp) const
{
	REQUIRE(NULL != fp);
	for (unsigned i = 0; i < 4; i++) {
		for (unsigned j = 0; j < 4; j++) {
			fprintf(fp, "%f ", _m[i][j]);
		}
		fprintf(fp, "\n");
	}
}

template<class T>
void tMatrix4<T>::scan(FILE *fp)
{
	REQUIRE(NULL != fp);
	for (unsigned i = 0; i < 4; i++) 
    {
		for (unsigned j = 0; j < 4; j++) 
        {
			value_type t;
			if (1 != treal<T>::scan(fp, &t)) 
            {
				GTB_ERROR("failed to scan matrix");
			}
			_m[i][j] = t;
		}
	}
}


template<class T>
ostream &operator<<(ostream &out, const tMatrix4<T> &m)
{
	for (unsigned i = 0; i < 4; i++) {
		for (unsigned j = 0; j < 4; j++) {
			out << m[i][j] << " ";
		}
		out << "\n";
	}
	return out;
}

template
ostream &operator<<(ostream &out, const tMatrix4<float> &m);
template
ostream &operator<<(ostream &out, const tMatrix4<double> &m);


template<>
void tMatrix4<float>::load()
{
    tMatrix4 MT = transpose();
    glLoadMatrixf(MT.as_array());
}

template<>
void tMatrix4<double>::load()
{
    tMatrix4 MT = transpose();
    glLoadMatrixd(MT.as_array());
}

template<>
void tMatrix4<float>::ogl_multiply()
{
    tMatrix4 MT = transpose();
    glMultMatrixf(MT.as_array());
}

template<>
void tMatrix4<double>::ogl_multiply()
{
    tMatrix4 MT = transpose();
    glMultMatrixd(MT.as_array());
}

/*--------- Tools ---------------*/
template<class T>
tMatrix4<T> TranslationMatrix(T Tx, T Ty, T Tz)
{
    tMatrix4<T> M = tMatrix4<T>::MATRIX4_IDENTITY;
    M[0][3] = Tx;
    M[1][3] = Ty;
    M[2][3] = Tz;

    return M;
}

template
tMatrix4<float> TranslationMatrix(float Tx, float Ty, float Tz);
template
tMatrix4<double> TranslationMatrix(double Tx, double Ty, double Tz);

template<class T>
tMatrix4<T> TranslationMatrix(const tPoint3<T>& Tv)
{
    return TranslationMatrix(Tv[0], Tv[1], Tv[2]);
}

template
tMatrix4<float> TranslationMatrix(const tPoint3<float>& Tv);
template
tMatrix4<double> TranslationMatrix(const tPoint3<double>& Tv);


template<class T>
tMatrix4<T> ScaleMatrix(T Sx, T Sy, T Sz)
{
    tMatrix4<T> M = tMatrix4<T>::MATRIX4_IDENTITY;
    M[0][0] = Sx;
    M[1][1] = Sy;
    M[2][2] = Sz;

    return M;
}
template
tMatrix4<float> ScaleMatrix(float Sx, float Sy, float Sz);
template
tMatrix4<double> ScaleMatrix(double Sx, double Sy, double Sz);

template<>
const tMatrix4<float> tMatrix4<float>::MATRIX4_ZERO(
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f);


template<>
const tMatrix4<float> tMatrix4<float>::MATRIX4_IDENTITY(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f);

template<>
const tMatrix4<double> tMatrix4<double>::MATRIX4_ZERO(0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0);


template<>
const tMatrix4<double> tMatrix4<double>::MATRIX4_IDENTITY(1.0, 0.0, 0.0, 0.0,
                                       0.0, 1.0, 0.0, 0.0,
                                       0.0, 0.0, 1.0, 0.0,
                                       0.0, 0.0, 0.0, 1.0);


template class tMatrix4<float>;
template class tMatrix4<double>;

GTB_END_NAMESPACE
