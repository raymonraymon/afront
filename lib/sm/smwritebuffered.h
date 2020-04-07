/*
===============================================================================

  FILE:  smwritebuffered.h
  
  CONTENTS:
  
    Writes a Streaming Mesh with a "little-cache" aware greedy reordering of
    triangles that is subject to a constraint of maximal delay of triangles.

  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    15 January 2005 -- radically simplified and improved
    16 December 2004 -- initial version created once all other crap was done
  
===============================================================================
*/
#ifndef SMWRITE_BUFFERED_H
#define SMWRITE_BUFFERED_H

#include "smwriter.h"

class SMwriteBuffered : public SMwriter
{
public:
  // smwriter interface function implementations

  void add_comment(const char* comment);

  void set_nverts(int nverts);
  void set_nfaces(int nfaces);
  void set_boundingbox(const float* bb_min_f, const float* bb_max_f);

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx);

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smwritebuffered functions

  bool open(SMwriter* smwriter, int max_delay=100);

  SMwriteBuffered();
  ~SMwriteBuffered();

private:
  SMwriter* smwriter;
  int max_delay;

  void write_triangle_delayed();
};

#endif
