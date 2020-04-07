
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



#ifndef _SFCURVE_H
#define _SFCURVE_H


class SpaceFillingCurve {
    public:
    SpaceFillingCurve(int x, int y, int z) {
		dim[0] = x;
		dim[1] = y;
		dim[2] = z;
    }

    virtual bool Done() const = 0;
    virtual void operator++() = 0;
    virtual void Get(int &x, int &y, int &z) = 0;

    protected:
    int dim[3];
};


class FlatCurve : public SpaceFillingCurve {
    public:
    FlatCurve(int x, int y, int z) : SpaceFillingCurve(x,y,z) {
		current[0] = 0;
		current[1] = 0;
		current[2] = 0;
    }

    bool Done() const {
		return (current[2] >= dim[2]);
    }


    void operator++() {

		// increment x
		current[0]++;

		// roll over into y
		current[1] += current[0] / dim[0];
		current[0]  = current[0] % dim[0];

		// roll over into z
		current[2] += current[1] / dim[1];
		current[1]  = current[1] % dim[1];
    }

    void Get(int &x, int &y, int &z) {
		x = current[0];
		y = current[1];
		z = current[2];
    }
    

    private:
    int current[3];
};

class ZCurve : public SpaceFillingCurve {
    public:

    ZCurve(int x, int y, int z) : SpaceFillingCurve(x,y,z) {
		index = 0;
		returned_indices = 0;
    }


    bool Done() const {
		return (returned_indices >= ((int64)dim[0]*(int64)dim[1]*(int64)dim[2]));
    }

    void operator++() {
		returned_indices++;
		index++;
    }

    void Get(int &x, int &y, int &z) {
		while (1) {
			i2xyz(index, x,y,z);
			if (x<dim[0] && y<dim[1] && z<dim[2])
				return;
			index++;
		}
    }


    static void xyz2i(int x, int y, int z, int64 &i) {

		const int64 c0 = 0x208208208208ULL;
		const int64 c1 = 0x041041041041ULL;
		const int64 c2 = 0x0c00c00c00c0ULL;
		const int64 c3 = 0x003003003003ULL;
		const int64 c4 = 0x00f00000f000ULL;
		const int64 c5 = 0x00000f00000fULL;
		const int64 c6 = 0x0000ff000000ULL;
		const int64 c7 = 0x0000000000ffULL;

		int64 xt, xt1, xt2;
		int64 yt, yt1, yt2;
		int64 zt, zt1, zt2;

		xt = x;                  yt = y;                  zt = z;

		xt1 = xt & c7;           yt1 = yt & c7;           zt1 = zt & c7;
		xt2 = (xt<<16) & c6;     yt2 = (yt<<16) & c6;     zt2 = (zt<<16) & c6;
		xt = xt1 | xt2;          yt = yt1 | yt2;          zt = zt1 | zt2; 

		xt1 = xt & c5;           yt1 = yt & c5;           zt1 = zt & c5;
		xt2 = (xt<<8) & c4;      yt2 = (yt<<8) & c4;      zt2 = (zt<<8) & c4;
		xt = xt1 | xt2;          yt = yt1 | yt2;          zt = zt1 | zt2; 

		xt1 = xt & c3;           yt1 = yt & c3;           zt1 = zt & c3;
		xt2 = (xt<<4) & c2;      yt2 = (yt<<4) & c2;      zt2 = (zt<<4) & c2;
		xt = xt1 | xt2;          yt = yt1 | yt2;          zt = zt1 | zt2; 

		xt1 = xt & c1;           yt1 = yt & c1;           zt1 = zt & c1;
		xt2 = (xt<<2) & c0;      yt2 = (yt<<2) & c0;      zt2 = (zt<<2) & c0;
		xt = xt1 | xt2;          yt = yt1 | yt2;          zt = zt1 | zt2; 


		i = (zt<<2) | (yt<<1) | xt;
    }


    static int64 x2i(int x) {

		const int64 c0 = 0x208208208208ULL;
		const int64 c1 = 0x041041041041ULL;
		const int64 c2 = 0x0c00c00c00c0ULL;
		const int64 c3 = 0x003003003003ULL;
		const int64 c4 = 0x00f00000f000ULL;
		const int64 c5 = 0x00000f00000fULL;
		const int64 c6 = 0x0000ff000000ULL;
		const int64 c7 = 0x0000000000ffULL;

		int64 xt, xt1, xt2;

		xt = x;                  

		xt1 = xt & c7;           
		xt2 = (xt<<16) & c6;     
		xt = xt1 | xt2;          

		xt1 = xt & c5;           
		xt2 = (xt<<8) & c4;      
		xt = xt1 | xt2;          

		xt1 = xt & c3;           
		xt2 = (xt<<4) & c2;      
		xt = xt1 | xt2;          

		xt1 = xt & c1;           
		xt2 = (xt<<2) & c0;      
		xt = xt1 | xt2;          

		return xt;
    }


    static void i2xyz(int64 i, int &x, int &y, int &z) {

		const int64 c0 = 0x208208208208ULL;
		const int64 c1 = 0x041041041041ULL;
		const int64 c2 = 0x0c00c00c00c0ULL;
		const int64 c3 = 0x003003003003ULL;
		const int64 c4 = 0x00f00000f000ULL;
		const int64 c5 = 0x00000f00000fULL;
		const int64 c6 = 0x0000ff000000ULL;
		const int64 c7 = 0x0000000000ffULL;

		int64 xt, xt1, xt2;
		int64 yt, yt1, yt2;
		int64 zt, zt1, zt2;

		xt = i>>0;               yt = i>>1;               zt = i>>2;

		xt1 = xt & c0;           yt1 = yt & c0;           zt1 = zt & c0;
		xt2 = xt & c1;           yt2 = yt & c1;           zt2 = zt & c1;
		xt  = (xt1>>2) | xt2;    yt  = (yt1>>2) | yt2;    zt  = (zt1>>2) | zt2;

		xt1 = xt & c2;           yt1 = yt & c2;           zt1 = zt & c2;
		xt2 = xt & c3;           yt2 = yt & c3;           zt2 = zt & c3;
		xt  = (xt1>>4) | xt2;    yt  = (yt1>>4) | yt2;    zt  = (zt1>>4) | zt2;

		xt1 = xt & c4;           yt1 = yt & c4;           zt1 = zt & c4;
		xt2 = xt & c5;           yt2 = yt & c5;           zt2 = zt & c5;
		xt  = (xt1>>8) | xt2;    yt  = (yt1>>8) | yt2;    zt  = (zt1>>8) | zt2;

		xt1 = xt & c6;           yt1 = yt & c6;           zt1 = zt & c6;
		xt2 = xt & c7;           yt2 = yt & c7;           zt2 = zt & c7;
		xt  = (xt1>>16) | xt2;   yt  = (yt1>>16) | yt2;   zt  = (zt1>>16) | zt2;

		x = (int)xt;
		y = (int)yt;
		z = (int)zt;
    }


    private:
    int64 index;
    int64 returned_indices;
};


#endif

