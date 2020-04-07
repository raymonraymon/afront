/*
===============================================================================

  FILE:  smreaddereferenced.h
  
  CONTENTS:
  
    Reads a Streaming Mesh and dereferences with the help of a hash to provide
    the t_pos_f field that contains pointers to the coordinates for the three
    vertices of the triangle. This makes it easy to read a streaming mesh as if
    it was "triangle soup".
    This filter can read streaming meshes in any order (e.g. post-order, pre-order,
    or in-order). Obviously the resulting streaming mesh will be in pre-order.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    25 July 2007 -- created while listening to a "Spoon" album on my old tivoli
  
===============================================================================
*/
#ifndef SMREAD_DEREFERENCED_H
#define SMREAD_DEREFERENCED_H

#include "smreader.h"

class SMreadDereferenced : public SMreader
{
public:
  // additional triangle variable
  float* t_pos_f[3];

  // smreader interface function implementations

  virtual SMevent read_element();
  virtual SMevent read_event();

  virtual void close(bool close_file=true); // virtual base function must be reimplemented

  // smreaddereferenced functions

  bool open(SMreader* smreader);

  SMreadDereferenced();
  virtual ~SMreadDereferenced();

private:
  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  SMreader* smreader;

  int read_input();
  int read_triangle();
};

#endif
