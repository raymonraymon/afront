/*
===============================================================================

  FILE:  smwriter_nil.h
  
  CONTENTS:
  
    A Streaming Mesh writer that is a dummy. It does not do anything but count.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    07 January 2006 -- created 
  
===============================================================================
*/
#ifndef SMWRITER_NIL_H
#define SMWRITER_NIL_H

#include "smwriter.h"

#include <stdio.h>

class SMwriter_nil : public SMwriter
{
public:

  // smwriter interface function implementations

  void add_comment(const char* comment) {};

  void set_boundingbox(const float* bb_min_f, const float* bb_max_f) {};

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx) {};

  void close(bool close_file=true, bool update_header=true) {v_count = -1; f_count = -1;};

  // smwriter_nil functions

  SMwriter_nil() {nverts = -1; nfaces = -1; v_count = 0; f_count = 0;};
  ~SMwriter_nil() {};
};

void SMwriter_nil::write_vertex(const float* v_pos_f)
{
  v_count++;
}
void SMwriter_nil::write_triangle(const int* t_idx, const bool* t_final)
{
  f_count++;
}
void SMwriter_nil::write_triangle(const int* t_idx)
{
  f_count++;
}
#endif
