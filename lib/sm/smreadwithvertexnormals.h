/*
===============================================================================

  FILE:  smreadwithvertexnormals.h
  
  CONTENTS:
  
    A filter for reading a Streaming Mesh with smooth vertex normals. It delays
    mesh elements until the three vertices of a triangle are finalized so that
    their smooth normals are available.

    The filter can deal with all stream orderings, pre-order, post-order, or
    in-order. The normal information is provided in the v_nor_f and the t_nor_f
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
  
    14 June 2007 -- created the day after i got official word of my LLNL offer 
  
===============================================================================
*/
#ifndef SMREAD_WITH_VERTEX_NORMALS_H
#define SMREAD_WITH_VERTEX_NORMALS_H

#include "smreaddereferenced.h"

#include <stdio.h>

class SMreadWithVertexNormals : public SMreadDereferenced
{
public:
  // additional vertex variable
  float v_nor_f[3];

  // additional triangle variable
  float* t_nor_f[3];

  // smreaddereferenced interface function implementations

  virtual SMevent read_element();
  virtual SMevent read_event();

  virtual void close(bool close_file=true); // virtual base function must be reimplemented

  // smreadwithvertexnormals functions

  bool open(SMreader* smreader);

  SMreadWithVertexNormals();
  ~SMreadWithVertexNormals();

private:
  SMreader* smreader;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  int read_input();
  int read_triangle();
};

#endif
