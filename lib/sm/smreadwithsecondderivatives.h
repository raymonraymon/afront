/*
===============================================================================

  FILE:  smreadwithsecondderivatives.h
  
  CONTENTS:
  
    A filter for reading a Streaming Mesh with smooth vertex normals alongsside
    first and second derivatives. It uses the SMreadWithVertexNormals filter and
    delays mesh elements until all the three vertices with smooth normals of a
    triangle are finalized so that their second derivatives are available.

    The filter can deal with all stream orderings, pre-order, post-order, or
    in-order. The gradient information is provided in the appropriately named
    field. The normal information is provided in the v_nor_f and the t_nor_f
    fields. The normals are not unit normals. You may call VecSelfNormalize3fv()
    on these normal fields.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    25 July 2007 -- created while unemployed and after the scary wine spill
  
===============================================================================
*/
#ifndef SMREAD_WITH_SECOND_DERIVATIVES_H
#define SMREAD_WITH_SECOND_DERIVATIVES_H

#include "smreadwithvertexnormals.h"

#include <stdio.h>

class SMreadWithSecondDerivatives : public SMreadWithVertexNormals
{
public:
  // additional vertex variable
  float v_dz_f[3];  // first derivatives v_dz_f[0] => dz/dx,  v_dz_f[1] = dz/dy,  v_dz_f[2] = 0 (unused)
  float v_d2z_f[3]; // second derivatives v_d2z_f[0] = d2z/d2x, v_d2z_f[1] = d2z/d2y, v_d2z_f[2] = d2z/dxdy

  // additional triangle variable
  float* t_dz_f[3]; 
  float* t_d2z_f[3];

  // smreader interface function implementations

  virtual SMevent read_element();
  virtual SMevent read_event();

  virtual void close(bool close_file=true); // virtual base function must be reimplemented

  // smreaderwithsecondderivatives functions

  bool open(SMreader* smreader);
  bool open(SMreadWithVertexNormals* smreadwithvertexnormals);

  SMreadWithSecondDerivatives();
  ~SMreadWithSecondDerivatives();

private:
  SMreadWithVertexNormals* smreadwithvertexnormals;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  int read_input();
  int read_triangle();
};

#endif
