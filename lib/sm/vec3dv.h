/*
===============================================================================

  FILE:  vec3dv.h
  
  CONTENTS:
  
    inlined functions for common vec3dv operations
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu  (with ideas stolen from Kenny Hoff)
  
  COPYRIGHT:
  
    copyright (C) 2000-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    11 January 2006 -- port from float-precision to double-precision
    12 August 2004 -- added copy with endian swap
    28 June 2000 -- created for GI 2000 submission
  
===============================================================================
*/
#ifndef VEC3DV_H
#define VEC3DV_H

#include <math.h>

inline void VecZero3dv(double v[3]);

inline void VecSet3dv(double v[3], double x, double y, double z);
inline void VecSet3dv_swap_endian(double v[3], double x, double y, double z);

inline void VecCopy3dv(double v[3], const double a[3]);
inline void VecCopy3dv_swap_endian(double v[3], const double a[3]);
inline void VecCopy3dv(double v[3], const float a[3]);

inline double VecMin3dv(const double v[3]);
inline double VecMax3dv(const double v[3]);

inline void VecScalarMult3dv(double v[3], const double a[3], double s);

inline void VecSelfScalarMult3dv(double v[3], double s);
inline void VecSelfScalarDiv3dv(double v[3], double s);
inline void VecSelfNegate3dv(double v[3]);

inline void VecSelfAdd3dv(double v[3], const double a[3]);
inline void VecSelfSubtract3dv(double v[3], const double a[3]);

inline void VecAdd3dv(double v[3], const double a[3], double b0, double b1, double b2);
inline void VecAdd3dv(double v[3], const double a[3], const double b[3]);
inline void VecAdd3dv(double v[3], const int a[3], const int b[3]);
inline void VecAdd3dv(double v[3], const double a[3], const double b[3], const double c[3]);
inline void VecAdd3dv(double v[3], const int a[3], const int b[3], const int c[3]);

inline void VecSubtract3dv(double v[3], const double a[3], double b0, double b1, double b2);
inline void VecSubtract3dv(double v[3], const double a[3], const double b[3]);
inline void VecSubtract3dv(double v[3], const int a[3], const int b[3]);

inline void VecAbsDiff3dv(double v[3], const double a[3], const double b[3]);
inline void VecScalarDiv3dv(double v[3], const double a[3], double s);

inline void VecParallelogram3dv(double v[3], const double a[3], const double b[3], const double c[3]);

inline void VecAddScalarMult3dv(double v[3], const double a[3], const double b[3], double s);
inline void VecSelfAddScalarMult3dv(double v[3], const double a[3], double s);

inline int VecEqual3dv(const double a[3], const double b[3]);

inline bool VecSelfNormalize3dv(double v[3]);
inline bool VecNormalize3dv(double v[3], const double a[3]);

inline double VecLength3dv(const double v[3]);
inline double VecDistance3dv(const double a[3], const double b[3]);
inline double VecSquaredLength3dv(const double v[3]);
inline double VecDotProd3dv(const double a[3], const double b[3]);

inline void VecCrossProd3dv(double v[3], const double a[3], const double b[3]);

inline void VecUpdateMinMax3dv(double min[3], double max[3], const double v[3]);
inline void VecUpdateMin3dv(double min[3], const double v[3]);
inline void VecUpdateMax3dv(double max[3], const double v[3]);
inline void VecUpdateMin3dv(double min[3], double x, double y, double z);
inline void VecUpdateMax3dv(double max[3], double x, double y, double z);

inline void VecClamp3dv(double v[3], const double min[3], const double max[3], const double a[3]);
inline void VecSelfClamp3dv(double v[3], const double min[3], const double max[3]);

inline void VecCcwNormal3dv(double n[3], const double a[3], const double b[3], const double c[3]);
inline void VecCwNormal3dv(double n[3], const double a[3], const double b[3], const double c[3]);

inline bool VecCcwNormNormal3dv(double n[3], const double a[3], const double b[3], const double c[3]);
inline bool VecCwNormNormal3dv(double n[3], const double a[3], const double b[3], const double c[3]);

inline bool VecCcwNormNormal3dv(double n[3], const int a[3], const int b[3], const int c[3]);
inline bool VecCwNormNormal3dv(double n[3], const int a[3], const int b[3], const int c[3]);

inline double VecCcwNormNormalReturnLength3dv(double n[3], const double a[3], const double b[3], const double c[3]);
inline double VecCwNormNormalReturnLength3dv(double n[3], const double a[3], const double b[3], const double c[3]);

inline double VecCcwNormNormalReturnLength3dv(double n[3], const int a[3], const int b[3], const int c[3]);
inline double VecCwNormNormalReturnLength3dv(double n[3], const int a[3], const int b[3], const int c[3]);

inline void VecCarthesianToSpherical3dv(double v[3], const double a[3]);
inline void VecSphericalToCarthesian3dv(double v[3], const double a[3]);

inline void VecZero3dv(double v[3])
{
  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = 0.0f;
}

inline void VecSet3dv(double v[3], double x, double y, double z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

inline void VecSet3dv_swap_endian(double v[3], double x, double y, double z)
{
  ((char*)v)[0] = ((const char*)&x)[7];
  ((char*)v)[1] = ((const char*)&x)[6];
  ((char*)v)[2] = ((const char*)&x)[5];
  ((char*)v)[3] = ((const char*)&x)[4];
  ((char*)v)[4] = ((const char*)&x)[3];
  ((char*)v)[5] = ((const char*)&x)[2];
  ((char*)v)[6] = ((const char*)&x)[1];
  ((char*)v)[7] = ((const char*)&x)[0];

  ((char*)v)[8] = ((const char*)&y)[7];
  ((char*)v)[9] = ((const char*)&y)[6];
  ((char*)v)[10] = ((const char*)&y)[5];
  ((char*)v)[11] = ((const char*)&y)[4];
  ((char*)v)[12] = ((const char*)&y)[3];
  ((char*)v)[13] = ((const char*)&y)[2];
  ((char*)v)[14] = ((const char*)&y)[1];
  ((char*)v)[15] = ((const char*)&y)[0];

  ((char*)v)[16] = ((const char*)&z)[7];
  ((char*)v)[17] = ((const char*)&z)[6];
  ((char*)v)[18] = ((const char*)&z)[5];
  ((char*)v)[19] = ((const char*)&z)[4];
  ((char*)v)[20] = ((const char*)&z)[3];
  ((char*)v)[21] = ((const char*)&z)[2];
  ((char*)v)[22] = ((const char*)&z)[1];
  ((char*)v)[23] = ((const char*)&z)[0];
}

inline void VecCopy3dv(double v[3], const double a[3])
{
  v[0] = a[0];
  v[1] = a[1];
  v[2] = a[2];
}

inline void VecCopy3dv(double v[3], const float a[3])
{
  v[0] = a[0];
  v[1] = a[1];
  v[2] = a[2];
}

inline void VecCopy3dv_swap_endian(double v[3], const double a[3])
{
  ((char*)v)[0] = ((const char*)a)[7];
  ((char*)v)[1] = ((const char*)a)[6];
  ((char*)v)[2] = ((const char*)a)[5];
  ((char*)v)[3] = ((const char*)a)[4];
  ((char*)v)[4] = ((const char*)a)[3];
  ((char*)v)[5] = ((const char*)a)[2];
  ((char*)v)[6] = ((const char*)a)[1];
  ((char*)v)[7] = ((const char*)a)[0];

  ((char*)v)[8] = ((const char*)a)[15];
  ((char*)v)[9] = ((const char*)a)[14];
  ((char*)v)[10] = ((const char*)a)[13];
  ((char*)v)[11] = ((const char*)a)[12];
  ((char*)v)[12] = ((const char*)a)[11];
  ((char*)v)[13] = ((const char*)a)[10];
  ((char*)v)[14] = ((const char*)a)[9];
  ((char*)v)[15] = ((const char*)a)[8];

  ((char*)v)[16] = ((const char*)a)[23];
  ((char*)v)[17] = ((const char*)a)[22];
  ((char*)v)[18] = ((const char*)a)[21];
  ((char*)v)[19] = ((const char*)a)[20];
  ((char*)v)[20] = ((const char*)a)[19];
  ((char*)v)[21] = ((const char*)a)[18];
  ((char*)v)[22] = ((const char*)a)[17];
  ((char*)v)[23] = ((const char*)a)[16];
}

inline double VecMin3dv(const double v[3])
{
  if (v[0] < v[1]) return ((v[0] < v[2]) ? v[0] : v[2]);
  else return ((v[1] < v[2]) ? v[1] : v[2]);
}

inline double VecMax3dv(const double v[3])
{
  if (v[0] > v[1]) return ((v[0] > v[2]) ? v[0] : v[2]);
  else return ((v[1] > v[2]) ? v[1] : v[2]);
}

inline void VecScalarMult3dv(double v[3], const double a[3], double s)
{
  v[0] = a[0]*s;
  v[1] = a[1]*s;
  v[2] = a[2]*s;
}

inline void VecSelfScalarMult3dv(double v[3], double s)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

inline void VecSelfScalarDiv3dv(double v[3], double s)
{
  v[0] /= s;
  v[1] /= s;
  v[2] /= s;
}

inline void VecSelfNegate3dv(double v[3])
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

inline void VecSelfAdd3dv(double v[3], const double a[3])
{
  v[0] += a[0];
  v[1] += a[1];
  v[2] += a[2];
}

inline void VecSelfSubtract3dv(double v[3], const double a[3])
{
  v[0] -= a[0];
  v[1] -= a[1];
  v[2] -= a[2];
}

inline void VecAdd3dv(double v[3], const double a[3], double b0, double b1, double b2)
{
  v[0] = a[0] + b0;
  v[1] = a[1] + b1;
  v[2] = a[2] + b2;
}

inline void VecAdd3dv(double v[3], const double a[3], const double b[3])
{
  v[0] = a[0] + b[0];
  v[1] = a[1] + b[1];
  v[2] = a[2] + b[2];
}

inline void VecAdd3dv(double v[3], const int a[3], const int b[3])
{
  v[0] = ((double)(a[0])) + ((double)(b[0]));
  v[1] = ((double)(a[1])) + ((double)(b[1]));
  v[2] = ((double)(a[2])) + ((double)(b[2]));
}

inline void VecAdd3dv(double v[3], const double a[3], const double b[3], const double c[3])
{
  v[0] = a[0] + b[0] + c[0];
  v[1] = a[1] + b[1] + c[1];
  v[2] = a[2] + b[2] + c[2];
}

inline void VecAdd3dv(double v[3], const int a[3], const int b[3], const int c[3])
{
  v[0] = ((double)(a[0])) + ((double)(b[0])) + ((double)(c[0]));
  v[1] = ((double)(a[1])) + ((double)(b[1])) + ((double)(c[1]));
  v[2] = ((double)(a[2])) + ((double)(b[2])) + ((double)(c[2]));
}

inline void VecSubtract3dv(double v[3], const double a[3], double b0, double b1, double b2)
{
  v[0] = a[0] - b0;
  v[1] = a[1] - b1;
  v[2] = a[2] - b2;
}

inline void VecSubtract3dv(double v[3], const double a[3], const double b[3])
{
  v[0] = a[0] - b[0];
  v[1] = a[1] - b[1];
  v[2] = a[2] - b[2];
}

inline void VecSubtract3dv(double v[3], const int a[3], const int b[3])
{
  v[0] = ((double)(a[0])) - ((double)(b[0]));
  v[1] = ((double)(a[1])) - ((double)(b[1]));
  v[2] = ((double)(a[2])) - ((double)(b[2]));
}

inline void VecAbsDiff3dv(double v[3], const double a[3], const double b[3])
{
  v[0] = a[0] - b[0];
  v[1] = a[1] - b[1];
  v[2] = a[2] - b[2];
  if (v[0] < 0.0f) v[0] = -v[0];
  if (v[1] < 0.0f) v[1] = -v[1];
  if (v[2] < 0.0f) v[2] = -v[2];
}

inline void VecScalarDiv3dv(double v[3], const double a[3], double s)
{
  v[0] = a[0]/s;
  v[1] = a[1]/s;
  v[2] = a[2]/s;
}

inline void VecParallelogram3dv(double v[3], const double a[3], const double b[3], const double c[3])
{
  v[0] = a[0] - b[0] + c[0];
  v[1] = a[1] - b[1] + c[1];
  v[2] = a[2] - b[2] + c[2];
}

inline void VecAddScalarMult3dv(double v[3], const double a[3], const double b[3], double s)
{
  v[0] = a[0] + (s * b[0]);
  v[1] = a[1] + (s * b[1]);
  v[2] = a[2] + (s * b[2]);
}

inline void VecSelfAddScalarMult3dv(double v[3], const double a[3], double s)
{
  v[0] += (s * a[0]);
  v[1] += (s * a[1]);
  v[2] += (s * a[2]);
}

inline double VecLength3dv(const double v[3])
{
  return( sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) );
}

inline double VecDistance3dv(const double a[3], const double b[3])
{
  return( sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2])) );
}

inline double VecSquaredLength3dv(const double v[3])
{
  return( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

inline bool VecNormalize3dv(double v[3], const double a[3])
{
  double l = VecLength3dv(a);
  if (l)
  {
    v[0] = a[0]/l;
    v[1] = a[1]/l;
    v[2] = a[2]/l;
    return true;
  }
  return false;
}

inline bool VecSelfNormalize3dv(double v[3])
{
  double l = VecLength3dv(v);
  if (l)
  {
    v[0] /= l;
    v[1] /= l;
    v[2] /= l;
    return true;
  }
  return false;
}

inline double VecDotProd3dv(const double a[3], const double b[3])
{
  return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

inline void VecCrossProd3dv(double v[3], const double a[3], const double b[3]) // c = a X b
{ 
  VecSet3dv(v, a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

inline void VecUpdateMinMax3dv(double min[3], double max[3], const double v[3])
{
  if (v[0]<min[0]) min[0]=v[0]; else if (v[0]>max[0]) max[0]=v[0];
  if (v[1]<min[1]) min[1]=v[1]; else if (v[1]>max[1]) max[1]=v[1];
  if (v[2]<min[2]) min[2]=v[2]; else if (v[2]>max[2]) max[2]=v[2];
}

inline void VecUpdateMin3dv(double min[3], const double v[3])
{
  if (v[0]<min[0]) min[0]=v[0];
  if (v[1]<min[1]) min[1]=v[1];
  if (v[2]<min[2]) min[2]=v[2];
}

inline void VecUpdateMax3dv(double max[3], const double v[3])
{
  if (v[0]>max[0]) max[0]=v[0];
  if (v[1]>max[1]) max[1]=v[1];
  if (v[2]>max[2]) max[2]=v[2];
}

inline void VecUpdateMin3dv(double min[3], double x, double y, double z)
{
  if (x<min[0]) min[0]=x;
  if (y<min[1]) min[1]=y;
  if (z<min[2]) min[2]=z;
}

inline void VecUpdateMax3dv(double max[3], double x, double y, double z)
{
  if (x>max[0]) max[0]=x;
  if (y>max[1]) max[1]=y;
  if (z>max[2]) max[2]=z;
}

inline void VecClamp3dv(double v[3], const double min[3], const double max[3], const double a[3])
{
  if (a[0]<min[0]) v[0] = min[0]; else if (a[0]>max[0]) v[0] = max[0]; else v[0] = a[0];
  if (a[1]<min[1]) v[1] = min[1]; else if (a[1]>max[1]) v[1] = max[1]; else v[1] = a[1];
  if (a[2]<min[2]) v[2] = min[2]; else if (a[2]>max[2]) v[2] = max[2]; else v[2] = a[2];
}

inline void VecSelfClamp3dv(double v[3], const double min[3], const double max[3])
{
  if (v[0]<min[0]) v[0] = min[0]; else if (v[0]>max[0]) v[0] = max[0];
  if (v[1]<min[1]) v[1] = min[1]; else if (v[1]>max[1]) v[1] = max[1];
  if (v[2]<min[2]) v[2] = min[2]; else if (v[2]>max[2]) v[2] = max[2];
}

inline int VecEqual3dv(const double a[3], const double b[3])
{
  return( a[0]==b[0] && a[1]==b[1] && a[2]==b[2] );
}

inline void VecCcwNormal3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ab,ac);
}

inline void VecCwNormal3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ac,ab);
}

inline bool VecCcwNormNormal3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ab,ac);
  return VecSelfNormalize3dv(n);
}

inline bool VecCwNormNormal3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ac,ab);
  return VecSelfNormalize3dv(n);
}

inline bool VecCcwNormNormal3dv(double n[3], const int a[3], const int b[3], const int c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ab,ac);
  return VecSelfNormalize3dv(n);
}

inline bool VecCwNormNormal3dv(double n[3], const int a[3], const int b[3], const int c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ac,ab);
  return VecSelfNormalize3dv(n);
}

inline double VecCcwNormNormalReturnLength3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ab,ac);
  double l = VecLength3dv(n);
  VecSelfScalarDiv3dv(n,l);
  return l;
}

inline double VecCwNormNormalReturnLength3dv(double n[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ac,ab);
  double l = VecLength3dv(n);
  VecSelfScalarDiv3dv(n,l);
  return l;
}

inline double VecCcwNormNormalReturnLength3dv(double n[3], const int a[3], const int b[3], const int c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ab,ac);
  double l = VecLength3dv(n);
  VecSelfScalarDiv3dv(n,l);
  return l;
}

inline double VecCwNormNormalReturnLength3dv(double n[3], const int a[3], const int b[3], const int c[3])
{
  double ab[3], ac[3];
  VecSubtract3dv(ab,b,a);
  VecSubtract3dv(ac,c,a);
  VecCrossProd3dv(n,ac,ab);
  double l = VecLength3dv(n);
  VecSelfScalarDiv3dv(n,l);
  return l;
}

inline void VecCarthesianToSpherical3dv(double v[3], const double a[3])
{
  double r2 = a[0]*a[0]+a[1]*a[1];
  v[1] = atan2(a[0], a[1]);
  v[2] = atan2(sqrt(r2), a[2]);
  v[0] = sqrt(r2+a[2]*a[2]);
}

inline void VecSphericalToCarthesian3dv(double v[3], const double a[3])
{
  v[0] = a[0]*sin(a[2])*sin(a[1]);
  v[1] = a[0]*sin(a[2])*cos(a[1]);
  v[2] = a[0]*cos(a[2]);
}

#endif
