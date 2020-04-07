
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



#ifndef _SPLINE_H
#define _SPLINE_H


template <typename T>
    class TrivariateSpline {
    public:

    TrivariateSpline() { }
	
    void SetCoefsCatmullRom(real_type tau=0.5) {
		CoefsFromMatrix(gtb::tMatrix4<double>(-tau, 2-tau, tau-2, tau,
											  2*tau, tau-3, 3-2*tau, -tau,
											  -tau, 0, tau, 0,
											  0, 1, 0, 0));
    }

    void SetCoefsBSpline() {
        CoefsFromMatrix(gtb::tMatrix4<double>(-1/6.0,  3/6.0, -3/6.0, 1/6.0,
											  3/6.0, -6/6.0,  3/6.0, 0/6.0,
											  -3/6.0,  0/6.0,  3/6.0, 0/6.0,
											  1/6.0,  4/6.0,  1/6.0, 0/6.0));
    }



	// use two different methods for computing spline intervals - slower but tighter bounds (less recursion)
	Interval Bound(T xaspect, T yaspect, T zaspect, 
				   const T nbrs[4][4][4], const T coefs[4][4][4],
				   const int partial[3], const Interval xspan_s[3],
				   int dx_sign=0, int dy_sign=0, int dz_sign=0, bool twomethod=false) const {

		Interval min_region[3] = {xspan_s[0], xspan_s[1], xspan_s[2] };
		Interval max_region[3] = {xspan_s[0], xspan_s[1], xspan_s[2] };


		int psigns[3] = { dx_sign, dy_sign, dz_sign };

		int axis_collapses = 0;

		for (int axis=0; axis<3; axis++) {
			if (psigns[axis]>0) {
				min_region[axis].max = min_region[axis].min;
				max_region[axis].min = max_region[axis].max;
				axis_collapses |= 1<<axis;
			}

			if (psigns[axis]<0) {
				min_region[axis].min = min_region[axis].max;
				max_region[axis].max = max_region[axis].min;
				axis_collapses |= 1<<axis;
			}
		}


		Interval ret;
		switch (axis_collapses) {
		case (1|2|4): { // all axis collapsed - eval at the corners
			real_type min_pt[3] = { min_region[0][0], min_region[1][0], min_region[2][0] };
			real_type max_pt[3] = { max_region[0][0], max_region[1][0], max_region[2][0] };
			
			real_type nmin = EvalPartial(xaspect, yaspect, zaspect, nbrs, min_pt, partial);
			real_type nmax = EvalPartial(xaspect, yaspect, zaspect, nbrs, max_pt, partial);

			ret = Interval(nmin, nmax);
			break;
		}
		
		case (1|2):
			{ // collapsed x and y
				real_type min_x = min_region[0][0];
				real_type max_x = max_region[0][0];

				real_type min_y = min_region[1][0];
				real_type max_y = max_region[1][0];

				Interval min_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, min_x, min_y, min_region[2], partial);
				Interval max_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, max_x, max_y, max_region[2], partial);

				// Interval min_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, min_x, min_y, min_region[2], partial);
				// Interval max_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, max_x, max_y, max_region[2], partial);

				ret = Interval(min_int.min, max_int.max);
				break;
			}
		
		case (1|4):
			{ // collapsed x and z
				real_type min_x = min_region[0][0];
				real_type max_x = max_region[0][0];

				real_type min_z = min_region[2][0];
				real_type max_z = max_region[2][0];
		
				Interval min_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, min_x, min_region[1], min_z, partial);
				Interval max_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, max_x, max_region[1], max_z, partial);

				// Interval min_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, min_x, min_region[1], min_z, partial);
				// Interval max_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, max_x, max_region[1], max_z, partial);

				ret = Interval(min_int.min, max_int.max);
				break;
			}
		
		case (2|4):
			{ // collapsed y and z
				real_type min_y = min_region[1][0];
				real_type max_y = max_region[1][0];
		
				real_type min_z = min_region[2][0];
				real_type max_z = max_region[2][0];

				Interval min_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, min_region[0], min_y, min_z, partial);
				Interval max_int = EvalPartial(xaspect, yaspect, zaspect, nbrs, max_region[0], max_y, max_z, partial);

				// Interval min_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, min_region[0], min_y, min_z, partial);
				// Interval max_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, max_region[0], max_y, max_z, partial);

				ret = Interval(min_int.min, max_int.max);
				break;
			}

		case (1):
		case (2):
		case (4):
			{
				// we collapsed just one dimension
				Interval min_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, min_region, partial);
				Interval max_int = EvalPartial2(xaspect, yaspect, zaspect, coefs, max_region, partial);


				if (twomethod) {
					min_int = Interval::Intersect(min_int, EvalPartial(xaspect, yaspect, zaspect, nbrs, min_region, partial));
					max_int = Interval::Intersect(max_int, EvalPartial(xaspect, yaspect, zaspect, nbrs, max_region, partial));
				}

				ret = Interval(min_int.min, max_int.max);
				break;
			}
		
		default:
			// default - just evaluate full region
			ret = EvalPartial2(xaspect, yaspect, zaspect, coefs, xspan_s, partial);

			if (twomethod) {
				ret = Interval::Intersect(ret, EvalPartial(xaspect, yaspect, zaspect, nbrs, xspan_s, partial));
			}

			break;
		}


		return ret;
	}




	template <typename T2>
    T2 Eval(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const T2 x[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };
	
		T2 zpre[4];
		T2 ypre[4];
		T2 xpre[4];
		for (int i=0; i<4; i++) {
			zpre[i] = cubic(matrix[i], xs[2]);
			ypre[i] = cubic(matrix[i], xs[1]);
			xpre[i] = cubic(matrix[i], xs[0]);
		}
 

		T2 accum = 0;
		for (int i=0; i<4; i++) {
			T2 &xcomp = xpre[i]; // cubic(matrix[i], xs[0]);

			for (int j=0; j<4; j++) {
				T2 &ycomp = ypre[j]; // cubic(matrix[j], xs[1]);

				for (int k=0; k<4; k++) {
					T2 &zcomp = zpre[k]; // cubic(matrix[k], xs[2]);

					accum += p[i][j][k] * xcomp * ycomp * zcomp;

				}
			}
		}

		return accum;
    }


	// accumulate the coefficients of the trivariate polynomial, then evaluate that directly
	template <typename T2>
	T2 Eval2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const T2 x[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };

		T2 xpow[4] = { cube(xs[0]), square(xs[0]), xs[0], 1 };
		T2 ypow[4] = { cube(xs[1]), square(xs[1]), xs[1], 1 };
		T2 zpow[4] = { cube(xs[2]), square(xs[2]), xs[2], 1 };


		T2 accum = 0;
		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {
					T2 pow = xpow[i] * ypow[j] * zpow[k];
					accum += coefs[i][j][k] * pow;
				}
			}
		}

		return accum;
    }


	template <typename T2>
	T2 EvalPartial(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const T2 x[3], const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };


		typedef double mat[4][4];
		const mat* allmats[4] = { &matrix, &dmatrix, &ddmatrix, &dddmatrix };
		const mat *xmat = allmats[partials[0]];
		const mat *ymat = allmats[partials[1]];
		const mat *zmat = allmats[partials[2]];

	
		T2 zpre[4];
		T2 ypre[4];
		T2 xpre[4];
		for (int i=0; i<4; i++) {
			xpre[i] = cubic((*xmat)[i], xs[0]);
			ypre[i] = cubic((*ymat)[i], xs[1]);
			zpre[i] = cubic((*zmat)[i], xs[2]);
		}
 

		T2 accum = 0;
		for (int i=0; i<4; i++) {
			T2 xcomp = xpre[i];

			for (int j=0; j<4; j++) {
				T2 ycomp = ypre[j];

				for (int k=0; k<4; k++) {
					T2 zcomp = zpre[k];

					accum += p[i][j][k] * xcomp * ycomp * zcomp;

				}
			}
		}

		return accum * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	Interval EvalPartial(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const Interval &x, T y, T z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;


		typedef double mat[4][4];
		const mat* allmats[4] = { &matrix, &dmatrix, &ddmatrix, &dddmatrix };
		const mat *xmat = allmats[partials[0]];
		const mat *ymat = allmats[partials[1]];
		const mat *zmat = allmats[partials[2]];

	
		T ypre[4];
		T zpre[4];
		for (int i=0; i<4; i++) {
			ypre[i] = cubic((*ymat)[i], y*yscale);
			zpre[i] = cubic((*zmat)[i], z*zscale);
		}



		T precoefs[4] = { 0,0,0,0 };
		for (int i=0; i<4; i++) {

			for (int j=0; j<4; j++) {
				T ycomp = ypre[j];

				for (int k=0; k<4; k++) {
					T zcomp = zpre[k];

					precoefs[i] += p[i][j][k] * ycomp * zcomp;
				}
			}
		}

		// result is now <xpre,precoef> - collapse it down to be a single cubic in x
		T ncubic[4];
		for (int j=0; j<4; j++) {
			ncubic[j] = 0;
			for (int i=0; i<4; i++) {
				ncubic[j] += (*xmat)[i][j] * precoefs[i];
			}
		}

		Interval res = cubic(ncubic, x*xscale);

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }


	Interval EvalPartial(T xaspect, T yaspect, T zaspect, const T p[4][4][4], T x, const Interval &y, T z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;


		typedef double mat[4][4];
		const mat* allmats[4] = { &matrix, &dmatrix, &ddmatrix, &dddmatrix };
		const mat *xmat = allmats[partials[0]];
		const mat *ymat = allmats[partials[1]];
		const mat *zmat = allmats[partials[2]];

	
		T xpre[4];
		T zpre[4];
		for (int i=0; i<4; i++) {
			xpre[i] = cubic((*xmat)[i], x*xscale);
			zpre[i] = cubic((*zmat)[i], z*zscale);
		}



		T precoefs[4] = { 0,0,0,0 };
		for (int i=0; i<4; i++) {
			T xcomp = xpre[i];

			for (int j=0; j<4; j++) {

				for (int k=0; k<4; k++) {
					T zcomp = zpre[k];

					precoefs[j] += p[i][j][k] * xcomp * zcomp;
				}
			}
		}

		// result is now <ypre,precoef> - collapse it down to be a single cubic in y
		T ncubic[4];
		for (int j=0; j<4; j++) {
			ncubic[j] = 0;
			for (int i=0; i<4; i++) {
				ncubic[j] += (*ymat)[i][j] * precoefs[i];
			}
		}

		Interval res = cubic(ncubic, y*yscale);

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	Interval EvalPartial(T xaspect, T yaspect, T zaspect, const T p[4][4][4], T x, T y, const Interval &z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		typedef double mat[4][4];
		const mat* allmats[4] = { &matrix, &dmatrix, &ddmatrix, &dddmatrix };

		const mat *xmat = allmats[partials[0]];
		const mat *ymat = allmats[partials[1]];
		const mat *zmat = allmats[partials[2]];

	
		T xpre[4];
		T ypre[4];
		for (int i=0; i<4; i++) {
			xpre[i] = cubic((*xmat)[i], x*xscale);
			ypre[i] = cubic((*ymat)[i], y*yscale);
		}



		T precoefs[4] = { 0,0,0,0 };
		for (int i=0; i<4; i++) {
			T xcomp = xpre[i];

			for (int j=0; j<4; j++) {
				T ycomp = ypre[j];

				for (int k=0; k<4; k++) {
					
					precoefs[k] += p[i][j][k] * xcomp * ycomp;
				}
			}
		}

		// result is now <zpre,precoef> - collapse it down to be a single cubic in z
		T ncubic[4];
		for (int j=0; j<4; j++) {
			ncubic[j] = 0;
			for (int i=0; i<4; i++) {
				ncubic[j] += (*zmat)[i][j] * precoefs[i];
			}
		}

		Interval res = cubic(ncubic, z*zscale);

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	Interval EvalPartial2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const Interval &x, T y, T z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T xs[3] = { 0, y*yscale, z*zscale };
		T pows[3][4];

		for (int i=1; i<3; i++) {
			switch (partials[i]) {
			case 0:
				pows[i][0] = cube(xs[i]);
				pows[i][1] = square(xs[i]);
				pows[i][2] = xs[i];
				pows[i][3] = 1;
				break;

			case 1:
				pows[i][0] = 3*square(xs[i]);
				pows[i][1] = 2*xs[i];
				pows[i][2] = 1;
				pows[i][3] = 0;
				break;

			case 2:
				pows[i][0] = 6*xs[i];
				pows[i][1] = 2;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;

			case 3:
				pows[i][0] = 6;
				pows[i][1] = 0;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;
			}
		}


 		T ncubic[4] = { 0,0,0,0 };
		for (int i=0; i<4-partials[0]; i++) {
			for (int j=0; j<4-partials[1]; j++) {
				for (int k=0; k<4-partials[2]; k++) {
					ncubic[i] += coefs[i][j][k] * pows[1][j] * pows[2][k];
				}
			}
		}


		// the result is now the partial[0]'th x derivative of ncubic(x)
		Interval res; 
		switch (partials[0]) {
		case 0:
			res = cubic(ncubic, x*xscale);
			break;

		case 1:
			ncubic[0] *= 3;
			ncubic[1] *= 2;
			res = quadratic(ncubic, x*xscale);
			break;

		case 2:
			ncubic[0] *= 6;
			ncubic[1] *= 2;
			res = linear(ncubic, x*xscale);
			break;

		case 3:
			res = ncubic[0]*6;
			break;
		}

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	Interval EvalPartial2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], T x, const Interval &y, T z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T xs[3] = { x*xscale, 0, z*zscale };
		T pows[3][4];

		for (int i=0; i<3; i+=2) {
			switch (partials[i]) {
			case 0:
				pows[i][0] = cube(xs[i]);
				pows[i][1] = square(xs[i]);
				pows[i][2] = xs[i];
				pows[i][3] = 1;
				break;

			case 1:
				pows[i][0] = 3*square(xs[i]);
				pows[i][1] = 2*xs[i];
				pows[i][2] = 1;
				pows[i][3] = 0;
				break;

			case 2:
				pows[i][0] = 6*xs[i];
				pows[i][1] = 2;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;

			case 3:
				pows[i][0] = 6;
				pows[i][1] = 0;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;
			}
		}


 		T ncubic[4] = { 0,0,0,0 };
		for (int i=0; i<4-partials[0]; i++) {
			for (int j=0; j<4-partials[1]; j++) {
				for (int k=0; k<4-partials[2]; k++) {
					ncubic[j] += coefs[i][j][k] * pows[0][i] * pows[2][k];
				}
			}
		}


		// the result is now the partial[0]'th y derivative of ncubic(y)
		Interval res; 
		switch (partials[0]) {
		case 0:
			res = cubic(ncubic, y*yscale);
			break;

		case 1:
			ncubic[0] *= 3;
			ncubic[1] *= 2;
			res = quadratic(ncubic, y*yscale);
			break;

		case 2:
			ncubic[0] *= 6;
			ncubic[1] *= 2;
			res = linear(ncubic, y*yscale);
			break;

		case 3:
			res = ncubic[0]*6;
			break;
		}

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	Interval EvalPartial2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], T x, T y, const Interval &z, const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T xs[3] = { x*xscale, y*yscale, 0 };
		T pows[3][4];

		for (int i=0; i<2; i++) {
			switch (partials[i]) {
			case 0:
				pows[i][0] = cube(xs[i]);
				pows[i][1] = square(xs[i]);
				pows[i][2] = xs[i];
				pows[i][3] = 1;
				break;

			case 1:
				pows[i][0] = 3*square(xs[i]);
				pows[i][1] = 2*xs[i];
				pows[i][2] = 1;
				pows[i][3] = 0;
				break;

			case 2:
				pows[i][0] = 6*xs[i];
				pows[i][1] = 2;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;

			case 3:
				pows[i][0] = 6;
				pows[i][1] = 0;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;
			}
		}


 		T ncubic[4] = { 0,0,0,0 };
		for (int i=0; i<4-partials[0]; i++) {
			for (int j=0; j<4-partials[1]; j++) {
				for (int k=0; k<4-partials[2]; k++) {
					ncubic[k] += coefs[i][j][k] * pows[0][i] * pows[1][j];
				}
			}
		}


		// the result is now the partial[0]'th z derivative of ncubic(z)
		Interval res; 
		switch (partials[0]) {
		case 0:
			res = cubic(ncubic, z*zscale);
			break;

		case 1:
			ncubic[0] *= 3;
			ncubic[1] *= 2;
			res = quadratic(ncubic, z*zscale);
			break;

		case 2:
			ncubic[0] *= 6;
			ncubic[1] *= 2;
			res = linear(ncubic, z*zscale);
			break;

		case 3:
			res = ncubic[0]*6;
			break;
		}

		return res * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }





	template <typename T2>
	T2 EvalPartial2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const T2 x[3], const int partials[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };

		T2 pows[3][4];

		for (int i=0; i<3; i++) {
			switch (partials[i]) {
			case 0:
				pows[i][0] = cube(xs[i]);
				pows[i][1] = square(xs[i]);
				pows[i][2] = xs[i];
				pows[i][3] = 1;
				break;

			case 1:
				pows[i][0] = 3*square(xs[i]);
				pows[i][1] = 2*xs[i];
				pows[i][2] = 1;
				pows[i][3] = 0;
				break;

			case 2:
				pows[i][0] = 6*xs[i];
				pows[i][1] = 2;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;

			case 3:
				pows[i][0] = 6;
				pows[i][1] = 0;
				pows[i][2] = 0;
				pows[i][3] = 0;
				break;
			}
		}


		T2 accum = 0;
		for (int i=0; i<4-partials[0]; i++) {
			for (int j=0; j<4-partials[1]; j++) {
				for (int k=0; k<4-partials[2]; k++) {
					accum += coefs[i][j][k] * pows[0][i] * pows[1][j] * pows[2][k];
				}
			}
		}

		return accum * pow(xscale, partials[0]) * pow(yscale, partials[1]) * pow(zscale, partials[2]);
    }



	template <typename T2>
    void Gradient(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const T2 x[3], T2 gradient[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };
	
		T2 zpre[4];
		T2 dzpre[4];
		T2 ypre[4];
		T2 dypre[4];
		T2 xpre[4];
		T2 dxpre[4];
		for (int i=0; i<4; i++) {
			zpre[i] = cubic(matrix[i], xs[2]);
			dzpre[i] = quadratic(&dmatrix[i][1], xs[2]);

			ypre[i] = cubic(matrix[i], xs[1]);
			dypre[i] = quadratic(&dmatrix[i][1], xs[1]);

			xpre[i] = cubic(matrix[i], xs[0]);
			dxpre[i] = quadratic(&dmatrix[i][1], xs[0]);
		}
 

		gradient[0] = 0;
		gradient[1] = 0;
		gradient[2] = 0;

		for (int i=0; i<4; i++) {
			T2 &xcomp = xpre[i]; // cubic(matrix[i], xs[0]);
			T2 &dxcomp = dxpre[i];

			for (int j=0; j<4; j++) {
				T2 &ycomp = ypre[j]; // cubic(matrix[j], xs[1]);
				T2 &dycomp = dypre[j];

				for (int k=0; k<4; k++) {
					T2 &zcomp = zpre[k]; // cubic(matrix[k], xs[2]);
					T2 &dzcomp = dzpre[k];

					gradient[0] += p[i][j][k] * dxcomp * ycomp * zcomp;
					gradient[1] += p[i][j][k] * xcomp * dycomp * zcomp;
					gradient[2] += p[i][j][k] * xcomp * ycomp * dzcomp;
				}
			}
		}

		gradient[0] *= xscale;
		gradient[1] *= yscale;
		gradient[2] *= zscale;
	}

	// accumulate the coefficients of the trivariate polynomial, then evaluate that directly
	template <typename T2>
    void Gradient2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const T2 x[3], T2 gradient[3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };

		T2 xpow[4] = { cube(xs[0]), square(xs[0]), xs[0], 1 };
		T2 dxpow[4] = { 3*square(xs[0]), 2*xs[0], 1, 0 };

		T2 ypow[4] = { cube(xs[1]), square(xs[1]), xs[1], 1 };
		T2 dypow[4] = { 3*square(xs[1]), 2*xs[1], 1, 0 };

		T2 zpow[4] = { cube(xs[2]), square(xs[2]), xs[2], 1 };
		T2 dzpow[4] = { 3*square(xs[2]), 2*xs[2], 1, 0 };


		gradient[0] = 0;
		gradient[1] = 0;
		gradient[2] = 0;
		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {

					gradient[0] += coefs[i][j][k] * dxpow[i] *  ypow[j] *  zpow[k];
					gradient[1] += coefs[i][j][k] *  xpow[i] * dypow[j] *  zpow[k];
					gradient[2] += coefs[i][j][k] *  xpow[i] *  ypow[j] * dzpow[k];
				}
			}
		}

		gradient[0] *= xscale;
		gradient[1] *= yscale;
		gradient[2] *= zscale;
    }



	template <typename T2>
    void Hessian(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const T2 x[3], gtb::tmat3<T2> &hessian) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };
	

		T2 zpre[4];
		T2 dzpre[4];
		T2 dzdzpre[4];
		T2 ypre[4];
		T2 dypre[4];
		T2 dydypre[4];
		T2 xpre[4];
		T2 dxpre[4];
		T2 dxdxpre[4];
		for (int i=0; i<4; i++) {
			zpre[i] = cubic(matrix[i], xs[2]);
			dzpre[i] = quadratic(&dmatrix[i][1], xs[2]);
			dzdzpre[i] = linear(&ddmatrix[i][2], xs[2]);

			ypre[i] = cubic(matrix[i], xs[1]);
			dypre[i] = quadratic(&dmatrix[i][1], xs[1]);
			dydypre[i] = linear(&ddmatrix[i][2], xs[1]);

			xpre[i] = cubic(matrix[i], xs[0]);
			dxpre[i] = quadratic(&dmatrix[i][1], xs[0]);
			dxdxpre[i] = linear(&ddmatrix[i][2], xs[0]);
		}
	
		hessian[0][0]=0;  hessian[0][1]=0;  hessian[0][2]=0;
		hessian[1][0]=0;  hessian[1][1]=0;  hessian[1][2]=0;
		hessian[2][0]=0;  hessian[2][1]=0;  hessian[2][2]=0;


		for (int i=0; i<4; i++) {
			T2 &xcomp = xpre[i]; // cubic(matrix[i], xs[0]);
			T2 &dxcomp = dxpre[i];
			T2 &dxdxcomp = dxdxpre[i];

			for (int j=0; j<4; j++) {
				T2 &ycomp = ypre[j]; // cubic(matrix[j], xs[1]);
				T2 &dycomp = dypre[j];
				T2 &dydycomp = dydypre[j];

				for (int k=0; k<4; k++) {
					T2 &zcomp = zpre[k]; // cubic(matrix[k], xs[2]);
					T2 &dzcomp = dzpre[k];
					T2 &dzdzcomp = dzdzpre[k];


					hessian[0][0] += p[i][j][k] * dxdxcomp * ycomp * zcomp;
					hessian[0][1] += p[i][j][k] * dxcomp * dycomp * zcomp;
					hessian[0][2] += p[i][j][k] * dxcomp * ycomp * dzcomp;

					hessian[1][1] += p[i][j][k] * xcomp * dydycomp * zcomp;
					hessian[1][2] += p[i][j][k] * xcomp * dycomp * dzcomp;

					hessian[2][2] += p[i][j][k] * xcomp * ycomp * dzdzcomp;
				}
			}
		}


		hessian[0][0] *= square(xscale);
		hessian[0][1] *= xscale * yscale;
		hessian[0][2] *= xscale * zscale;

		hessian[1][1] *= square(yscale);
		hessian[1][2] *= yscale * zscale;
		
		hessian[2][2] *= square(zscale);


		hessian[1][0] = hessian[0][1];
		hessian[2][0] = hessian[0][2];
		hessian[2][1] = hessian[1][2];

    }


	// accumulate the coefficients of the trivariate polynomial, then evaluate that directly
	template <typename T2>
    void Hessian2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const T2 x[3], gtb::tmat3<T2> &hessian) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };


		T2 xpow[4] = { cube(xs[0]), square(xs[0]), xs[0], 1 };
		T2 dxpow[4] = { 3*square(xs[0]), 2*xs[0], 1, 0 };
		T2 dxdxpow[4] = { 6*xs[0], 2, 0, 0 };

		T2 ypow[4] = { cube(xs[1]), square(xs[1]), xs[1], 1 };
		T2 dypow[4] = { 3*square(xs[1]), 2*xs[1], 1, 0 };
		T2 dydypow[4] = { 6*xs[1], 2, 0, 0 };

		T2 zpow[4] = { cube(xs[2]), square(xs[2]), xs[2], 1 };
		T2 dzpow[4] = { 3*square(xs[2]), 2*xs[2], 1, 0 };
		T2 dzdzpow[4] = { 6*xs[2], 2, 0, 0 };



		hessian[0][0]=0;  hessian[0][1]=0;  hessian[0][2]=0;
		hessian[1][0]=0;  hessian[1][1]=0;  hessian[1][2]=0;
		hessian[2][0]=0;  hessian[2][1]=0;  hessian[2][2]=0;

		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {

					hessian[0][0] += coefs[i][j][k] * dxdxpow[i] *    ypow[j] *    zpow[k];
					hessian[0][1] += coefs[i][j][k] *   dxpow[i] *   dypow[j] *    zpow[k];
					hessian[0][2] += coefs[i][j][k] *   dxpow[i] *    ypow[j] *   dzpow[k];

					hessian[1][1] += coefs[i][j][k] *    xpow[i] * dydypow[j] *    zpow[k];
					hessian[1][2] += coefs[i][j][k] *    xpow[i] *   dypow[j] *   dzpow[k];

					hessian[2][2] += coefs[i][j][k] *    xpow[i] *    ypow[j] * dzdzpow[k];
				}
			}
		}

		hessian[0][0] *= square(xscale);
		hessian[0][1] *= xscale * yscale;
		hessian[0][2] *= xscale * zscale;

		hessian[1][1] *= square(yscale);
		hessian[1][2] *= yscale * zscale;
		
		hessian[2][2] *= square(zscale);


		hessian[1][0] = hessian[0][1];
		hessian[2][0] = hessian[0][2];
		hessian[2][1] = hessian[1][2];
    }



	template <typename T2>
    void ThirdTensor(T xaspect, T yaspect, T zaspect, const T p[4][4][4], const T2 x[3], T2 third[3][3][3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };
	

		T2 zpre[4];
		T2 dzpre[4];
		T2 dzdzpre[4];
		T2 dzdzdzpre[4];
		T2 ypre[4];
		T2 dypre[4];
		T2 dydypre[4];
		T2 dydydypre[4];
		T2 xpre[4];
		T2 dxpre[4];
		T2 dxdxpre[4];
		T2 dxdxdxpre[4];
		for (int i=0; i<4; i++) {
			zpre[i] = cubic(matrix[i], xs[2]);
			dzpre[i] = quadratic(&dmatrix[i][1], xs[2]);
			dzdzpre[i] = linear(&ddmatrix[i][2], xs[2]);
			dzdzdzpre[i] = dddmatrix[i][3];

			ypre[i] = cubic(matrix[i], xs[1]);
			dypre[i] = quadratic(&dmatrix[i][1], xs[1]);
			dydypre[i] = linear(&ddmatrix[i][2], xs[1]);
			dydydypre[i] = dddmatrix[i][3];

			xpre[i] = cubic(matrix[i], xs[0]);
			dxpre[i] = quadratic(&dmatrix[i][1], xs[0]);
			dxdxpre[i] = linear(&ddmatrix[i][2], xs[0]);
			dxdxdxpre[i] = dddmatrix[i][3];
		}

		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					third[i][j][k] = 0;
				}
			}
		}

		for (int i=0; i<4; i++) {
			T2 &xcomp = xpre[i]; // cubic(matrix[i], xs[0]);
			T2 &dxcomp = dxpre[i];
			T2 &dxdxcomp = dxdxpre[i];
			T2 &dxdxdxcomp = dxdxdxpre[i];

			for (int j=0; j<4; j++) {
				T2 &ycomp = ypre[j]; // cubic(matrix[j], xs[1]);
				T2 &dycomp = dypre[j];
				T2 &dydycomp = dydypre[j];
				T2 &dydydycomp = dydydypre[j];

				for (int k=0; k<4; k++) {
					T2 &zcomp = zpre[k]; // cubic(matrix[k], xs[2]);
					T2 &dzcomp = dzpre[k];
					T2 &dzdzcomp = dzdzpre[k];
					T2 &dzdzdzcomp = dzdzdzpre[k];


					third[0][0][0] += p[i][j][k] * dxdxdxcomp * ycomp * zcomp;
					third[0][0][1] += p[i][j][k] * dxdxcomp * dycomp * zcomp;
					third[0][0][2] += p[i][j][k] * dxdxcomp * ycomp * dzcomp;

					third[0][1][1] += p[i][j][k] * dxcomp * dydycomp * zcomp;
					third[0][1][2] += p[i][j][k] * dxcomp * dycomp * dzcomp;

					third[0][2][2] += p[i][j][k] * dxcomp * ycomp * dzdzcomp;


					third[1][1][1] += p[i][j][k] * xcomp * dydydycomp * zcomp;
					third[1][1][2] += p[i][j][k] * xcomp * dydycomp * dzcomp;

					third[1][2][2] += p[i][j][k] * xcomp * dycomp * dzdzcomp;

					third[2][2][2] += p[i][j][k] * xcomp * ycomp * dzdzdzcomp;
				}
			}
		}

		third[0][0][0] *= cube(xscale);
		third[0][0][1] *= square(xscale) * yscale;
		third[0][0][2] *= square(xscale) * zscale;

		third[0][1][1] *= xscale * square(yscale);
		third[0][1][2] *= xscale * yscale * zscale;


		third[0][2][2] *= xscale * square(zscale);


		third[1][1][1] *= cube(yscale);
		third[1][1][2] *= square(yscale) * zscale;

		third[1][2][2] *= yscale * square(zscale);

		third[2][2][2] *= cube(zscale);


		// fill in the symmetry
		third[0][1][0] = third[0][0][1];
		third[0][2][0] = third[0][0][2];
		third[0][2][1] = third[0][1][2];

		third[1][0][0] = third[0][0][1];
		third[1][0][1] = third[0][1][1];
		third[1][0][2] = third[0][1][2];
		third[1][1][0] = third[0][1][1];
		third[1][2][0] = third[0][1][2];
		third[1][2][1] = third[1][1][2];


		third[2][0][0] = third[0][0][2];
		third[2][0][1] = third[0][1][2];
		third[2][0][2] = third[0][2][2];
		third[2][1][0] = third[0][1][2];
		third[2][1][1] = third[1][1][2];
		third[2][1][2] = third[1][2][2];
		third[2][2][0] = third[0][2][2];
		third[2][2][1] = third[1][2][2];

    }


	template <typename T2>
    void ThirdTensor2(T xaspect, T yaspect, T zaspect, const T coefs[4][4][4], const T2 x[3], T2 third[3][3][3]) const {

		T xscale = 1.0 / xaspect;
		T yscale = 1.0 / yaspect;
		T zscale = 1.0 / zaspect;

		T2 xs[3] = { x[0]*xscale, x[1]*yscale, x[2]*zscale };


		T2 xpow[4] = { cube(xs[0]), square(xs[0]), xs[0], 1 };
		T2 dxpow[4] = { 3*square(xs[0]), 2*xs[0], 1, 0 };
		T2 dxdxpow[4] = { 6*xs[0], 2, 0, 0 };
		T2 dxdxdxpow[4] = { 6, 0, 0, 0 };

		T2 ypow[4] = { cube(xs[1]), square(xs[1]), xs[1], 1 };
		T2 dypow[4] = { 3*square(xs[1]), 2*xs[1], 1, 0 };
		T2 dydypow[4] = { 6*xs[1], 2, 0, 0 };
		T2 dydydypow[4] = { 6, 0, 0, 0 };

		T2 zpow[4] = { cube(xs[2]), square(xs[2]), xs[2], 1 };
		T2 dzpow[4] = { 3*square(xs[2]), 2*xs[2], 1, 0 };
		T2 dzdzpow[4] = { 6*xs[2], 2, 0, 0 };
		T2 dzdzdzpow[4] = { 6, 0, 0, 0 };



		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				for (int k=0; k<3; k++) {
					third[i][j][k] = 0;
				}
			}
		}




		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {

					third[0][0][0] += coefs[i][j][k] * dxdxdxpow[i] * ypow[j] * zpow[k];
					third[0][0][1] += coefs[i][j][k] * dxdxpow[i] * dypow[j] * zpow[k];
					third[0][0][2] += coefs[i][j][k] * dxdxpow[i] * ypow[j] * dzpow[k];

					third[0][1][1] += coefs[i][j][k] * dxpow[i] * dydypow[j] * zpow[k];
					third[0][1][2] += coefs[i][j][k] * dxpow[i] * dypow[j] * dzpow[k];

					third[0][2][2] += coefs[i][j][k] * dxpow[i] * ypow[j] * dzdzpow[k];


					third[1][1][1] += coefs[i][j][k] * xpow[i] * dydydypow[j] * zpow[k];
					third[1][1][2] += coefs[i][j][k] * xpow[i] * dydypow[j] * dzpow[k];

					third[1][2][2] += coefs[i][j][k] * xpow[i] * dypow[j] * dzdzpow[k];

					third[2][2][2] += coefs[i][j][k] * xpow[i] * ypow[j] * dzdzdzpow[k];
				}
			}
		}


		third[0][0][0] *= cube(xscale);
		third[0][0][1] *= square(xscale) * yscale;
		third[0][0][2] *= square(xscale) * zscale;

		third[0][1][1] *= xscale * square(yscale);
		third[0][1][2] *= xscale * yscale * zscale;


		third[0][2][2] *= xscale * square(zscale);


		third[1][1][1] *= cube(yscale);
		third[1][1][2] *= square(yscale) * zscale;

		third[1][2][2] *= yscale * square(zscale);

		third[2][2][2] *= cube(zscale);


		// fill in the symmetry
		third[0][1][0] = third[0][0][1];
		third[0][2][0] = third[0][0][2];
		third[0][2][1] = third[0][1][2];

		third[1][0][0] = third[0][0][1];
		third[1][0][1] = third[0][1][1];
		third[1][0][2] = third[0][1][2];
		third[1][1][0] = third[0][1][1];
		third[1][2][0] = third[0][1][2];
		third[1][2][1] = third[1][1][2];


		third[2][0][0] = third[0][0][2];
		third[2][0][1] = third[0][1][2];
		third[2][0][2] = third[0][2][2];
		third[2][1][0] = third[0][1][2];
		third[2][1][1] = third[1][1][2];
		third[2][1][2] = third[1][2][2];
		third[2][2][0] = third[0][2][2];
		third[2][2][1] = third[1][2][2];

    }



    void CoefsFromMatrix(const gtb::tMatrix4<T> &m) {
		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				matrix[i][j] = m[j][i];
			}

			dmatrix[i][0] = 0;
			dmatrix[i][1] = 3*matrix[i][0];
			dmatrix[i][2] = 2*matrix[i][1];
			dmatrix[i][3] = 1*matrix[i][2];

			ddmatrix[i][0] = 0;
			ddmatrix[i][1] = 0;
			ddmatrix[i][2] = 6*matrix[i][0];
			ddmatrix[i][3] = 2*matrix[i][1];

			dddmatrix[i][0] = 0;
			dddmatrix[i][1] = 0;
			dddmatrix[i][2] = 0;
			dddmatrix[i][3] = 6*matrix[i][0];
		}
    }



	// accumulate the coefficients of the trivariate polynomial, then evaluate that directly
	void Precompute(const T p[4][4][4], T coefs[4][4][4]) const {

		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {
					coefs[i][j][k] = 0;
				}
			}
		}

		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {

					// multiply out the tensor product polynomail and accumulate the coefficients of the trivariate
					for (int i2=0; i2<4; i2++) {
						for (int j2=0; j2<4; j2++) {
							for (int k2=0; k2<4; k2++) {
								coefs[i2][j2][k2] += p[i][j][k] 
									* matrix[i][i2]
									* matrix[j][j2]
									* matrix[k][k2];
							}
						}
					}

				}
			}
		}
    }




    private:
    double matrix[4][4];
	double dmatrix[4][4];
	double ddmatrix[4][4];
	double dddmatrix[4][4];
};


#endif


