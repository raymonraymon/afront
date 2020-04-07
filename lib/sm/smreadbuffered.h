/*
===============================================================================

  FILE:  smreadbuffered.h
  
  CONTENTS:
  
    Reads a Streaming Mesh with a "little-cache" aware greedy reordering of
    triangles that is subject to a constraint of maximal delay of triangles.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    27 June 2005 -- created to provide functionality of SMwriteBuffered on read
  
===============================================================================
*/
#ifndef SMREAD_BUFFERED_H
#define SMREAD_BUFFERED_H

#include "smreader.h"

class SMreadBuffered : public SMreader
{
public:
  // smreader interface function implementations

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreadbuffered functions

  bool open(SMreader* smreader, int max_delay=100);

  SMreadBuffered();
  ~SMreadBuffered();

  void* wait; // points to the waiting queue of triangles (this is a hack for the sm_viewer)

private:
  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  SMreader* smreader;
  int max_delay;

  int read_triangle_delayed();
  void fill_triangle_buffer();
};

#endif
