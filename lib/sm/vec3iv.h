/*
===============================================================================

  FILE:  vec3iv.h
  
  CONTENTS:
  
    inlined functions for common vec3iv operations
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu  (with ideas stolen from Kenny Hoff)
  
  COPYRIGHT:
  
    copyright (C) 2000-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    12 August 2004 -- added copy with endian swap
    11 January 2003 -- fixed up for SIGGRAPH 2003
    28 June 2000 -- created for GI 2000 submission
  
===============================================================================
*/
#ifndef VEC3IV_H
#define VEC3IV_H

#include <math.h>

inline void VecZero3iv(int v[3]);

inline void VecSet3iv(int v[3], int x, int y, int z);
inline void VecSet3iv_swap_endian(int v[3], int x, int y, int z);

inline void VecCopy3iv(int v[3], const int a[3]);
inline void VecCopy3iv_swap_endian(int v[3], const int a[3]);

inline int VecMin3iv(const int v[3]);
inline int VecMax3iv(const int v[3]);

inline void VecSelfScalarMult3iv(int v[3], int s);
inline void VecSelfScalarDiv3iv(int v[3], int s);
inline void VecSelfNegate3iv(int v[3]);

inline void VecSelfAdd3iv(int v[3], const int a[3]);
inline void VecSelfSubtract3iv(int v[3], const int a[3]);

inline void VecSelfAdd3iv(int v[3], int s);
inline void VecSelfSubtract3iv(int v[3], int s);

inline void VecAdd3iv(int v[3], const int a[3], const int b[3]);
inline void VecAdd3iv(int v[3], const int a[3], const int b[3], const int c[3]);
inline void VecSubtract3iv(int v[3], const int a[3], const int b[3]);

inline int VecEqual3iv(const int a[3], const int b[3]);

inline void VecNormalize3iv(float v[3], const int a[3]);

inline float VecLength3iv(const int v[3]);
inline float VecLength3iv(int v0, int v1, int v2);
inline float VecDistance3fiv(const int a[3], const int b[3]);
inline int  VecSquaredLength3iv(const int v[3]);
inline int  VecDotProd3iv(const int a[3], const int b[3]);

inline void VecCrossProd3iv(int v[3], const int a[3], const int b[3]);

inline void VecUpdateMinMax3iv(int min[3], int max[3], const int v[3]);

inline void VecCcwNormal3iv(int n[3], const int a[3], const int b[3], const int c[3]);
inline void VecCwNormal3iv(int n[3], const int a[3], const int b[3], const int c[3]);

inline int VecCcwNormNormal3iv(float n[3], const int a[3], const int b[3], const int c[3]);
inline int VecCwNormNormal3iv(float n[3], const int a[3], const int b[3], const int c[3]);

inline void VecZero3iv(int v[3])
{
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
}

inline void VecSet3iv(int v[3], int x, int y, int z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void VecSet3iv_swap_endian(int v[3], int x, int y, int z)
{
  ((char*)v)[0] = ((const char*)&x)[3];
  ((char*)v)[1] = ((const char*)&x)[2];
  ((char*)v)[2] = ((const char*)&x)[1];
  ((char*)v)[3] = ((const char*)&x)[0];
  ((char*)v)[4] = ((const char*)&y)[3];
  ((char*)v)[5] = ((const char*)&y)[2];
  ((char*)v)[6] = ((const char*)&y)[1];
  ((char*)v)[7] = ((const char*)&y)[0];
  ((char*)v)[8] = ((const char*)&z)[3];
  ((char*)v)[9] = ((const char*)&z)[2];
  ((char*)v)[10] = ((const char*)&z)[1];
  ((char*)v)[11] = ((const char*)&z)[0];
}

inline void VecCopy3iv(int v[3], const int a[3])
{
  v[0] = a[0];
  v[1] = a[1];
  v[2] = a[2];
}

inline void VecCopy3iv_swap_endian(int v[3], const int a[3])
{
  ((char*)v)[0] = ((const char*)a)[3];
  ((char*)v)[1] = ((const char*)a)[2];
  ((char*)v)[2] = ((const char*)a)[1];
  ((char*)v)[3] = ((const char*)a)[0];
  ((char*)v)[4] = ((const char*)a)[7];
  ((char*)v)[5] = ((const char*)a)[6];
  ((char*)v)[6] = ((const char*)a)[5];
  ((char*)v)[7] = ((const char*)a)[4];
  ((char*)v)[8] = ((const char*)a)[11];
  ((char*)v)[9] = ((const char*)a)[10];
  ((char*)v)[10] = ((const char*)a)[9];
  ((char*)v)[11] = ((const char*)a)[8];
}

inline int VecMin3iv(const int v[3])
{
  if (v[0] < v[1]) return ((v[0] < v[2]) ? v[0] : v[2]);
  else return ((v[1] < v[2]) ? v[1] : v[2]);
}

inline int VecMax3iv(const int v[3])
{
  if (v[0] > v[1]) return ((v[0] > v[2]) ? v[0] : v[2]);
  else return ((v[1] > v[2]) ? v[1] : v[2]);
}

inline void VecSelfScalarMult3iv(int v[3], int s)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

inline void VecSelfScalarDiv3iv(int v[3], int s)
{
  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
}

inline void VecSelfNegate3iv(int v[3])
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

inline void VecSelfAdd3iv(int v[3], const int a[3])
{
  v[0] += a[0];
  v[1] += a[1];
  v[2] += a[2];
}

inline void VecSelfSubtract3iv(int v[3], const int a[3])
{
  v[0] -= a[0];
  v[1] -= a[1];
  v[2] -= a[2];
}

inline void VecSelfAdd3iv(int v[3], int s)
{
  v[0] += s;
  v[1] += s;
  v[2] += s;
}

inline void VecSelfSubtract3iv(int v[3], int s)
{
  v[0] -= s;
  v[1] -= s;
  v[2] -= s;
}

inline void VecAdd3iv(int v[3], const int a[3], const int b[3])
{
  v[0] = a[0] + b[0];
  v[1] = a[1] + b[1];
  v[2] = a[2] + b[2];
}

inline void VecAdd3iv(int v[3], const int a[3], const int b[3], const int c[3])
{
  v[0] = a[0] + b[0] + c[0];
  v[1] = a[1] + b[1] + c[1];
  v[2] = a[2] + b[2] + c[2];
}

inline void VecSubtract3iv(int v[3], const int a[3], const int b[3])
{
  v[0] = a[0] - b[0];
  v[1] = a[1] - b[1];
  v[2] = a[2] - b[2];
}

inline float VecLength3iv(const int v[3])
{
  return( (float)sqrt( (float)v[0]*(float)v[0] + (float)v[1]*(float)v[1] + (float)v[2]*(float)v[2]) );
}

inline float VecLength3iv(int v0, int v1, int v2)
{
  return( (float)sqrt( (float)v0*(float)v0 + (float)v1*(float)v1 + (float)v2*(float)v2) );
}

inline int VecSquaredLength3iv(const int v[3])
{
  return( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

inline float VecDistance3iv(const int a[3], const int b[3])
{
  return( (float)sqrt(((float)(a[0]-b[0]))*((float)(a[0]-b[0])) + ((float)(a[1]-b[1]))*((float)(a[1]-b[1])) + ((float)(a[2]-b[2]))*((float)(a[2]-b[2]))) );
}

inline void VecNormalize3iv(float v[3], const int a[3])
{
  float l = VecLength3iv(a);
  v[0] = float(a[0])/l;
  v[1] = float(a[1])/l;
  v[2] = float(a[2])/l;
}

inline int VecDotProd3iv(const int a[3], const int b[3])
{
  return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

inline void VecCrossProd3iv(int v[3], const int a[3], const int b[3]) // v = a X b
{ 
  VecSet3iv(v, a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

inline void VecUpdateMinMax3iv(int min[3], int max[3], const int v[3])
{
  if (v[0]<min[0]) min[0]=v[0]; else if (v[0]>max[0]) max[0]=v[0];
  if (v[1]<min[1]) min[1]=v[1]; else if (v[1]>max[1]) max[1]=v[1];
  if (v[2]<min[2]) min[2]=v[2]; else if (v[2]>max[2]) max[2]=v[2];
}

inline int VecEqual3iv(const int a[3], const int b[3])
{
  return( a[0]==b[0] && a[1]==b[1] && a[2]==b[2] );
}

inline void VecCcwNormal3iv(int n[3], const int a[3], const int b[3], const int c[3])
{
  int ab[3], ac[3];
  VecSubtract3iv(ab,b,a);
  VecSubtract3iv(ac,c,a);
  VecCrossProd3iv(n,ab,ac);
}

inline void VecCwNormal3iv(int n[3], const int a[3], const int b[3], const int c[3])
{
  int ab[3], ac[3];
  VecSubtract3iv(ab,b,a);
  VecSubtract3iv(ac,c,a);
  VecCrossProd3iv(n,ac,ab);
}

inline int VecCcwNormNormal3iv(float n[3], const int a[3], const int b[3], const int c[3])
{
  float l;
  int ab[3], ac[3], abac[3];
  VecSubtract3iv(ab,b,a);
  VecSubtract3iv(ac,c,a);
  VecCrossProd3iv(abac,ab,ac);
  if (abac[0] || abac[1] || abac[2])
  {
    while (abac[0] > 16384 || abac[0] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    while (abac[1] > 16384 || abac[1] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    while (abac[2] > 16384 || abac[2] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    l = (float)sqrt((double)VecSquaredLength3iv(abac));
    n[0] = float(abac[0])/l;
    n[1] = float(abac[1])/l;
    n[2] = float(abac[2])/l;
    return 1;
  }
  else
  {
    return 0;
  }
}

inline int VecCwNormNormal3iv(float n[3], const int a[3], const int b[3], const int c[3])
{
  float l;
  int ab[3], ac[3], abac[3];
  VecSubtract3iv(ab,b,a);
  VecSubtract3iv(ac,c,a);
  VecCrossProd3iv(abac,ac,ab);
  if (abac[0] || abac[1] || abac[2])
  {
    while (abac[0] > 16384 || abac[0] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    while (abac[1] > 16384 || abac[1] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    while (abac[2] > 16384 || abac[2] < -16384) 
    {
      abac[0] >>= 6;
      abac[1] >>= 6;
      abac[2] >>= 6;
    }
    l = (float)sqrt((double)VecSquaredLength3iv(abac));
    n[0] = float(abac[0])/l;
    n[1] = float(abac[1])/l;
    n[2] = float(abac[2])/l;
    return 1;
  }
  else
  {
    return 0;
  }
}

#endif
