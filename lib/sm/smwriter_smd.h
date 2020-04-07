/*
===============================================================================

  FILE:  smwriter_smd.h
  
  CONTENTS:
  
    Writes a compressed Streaming Mesh using a more involved compression scheme
    that utilized a small buffer from which it chooses its triangles. There is
    a user-defined constraint on the maximum that a triangle can be *delayed* in
    the buffer.
    
    There are two ways to specify this constraint. Either as an abolute number
    delay (for 10000 triangles set delay=10000) or as an adaptive delay that is a
    multiple of the current width (for 2 times the width set delay=-2).
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 May 2005 -- fixed a Microsoft bug (floating-point in Release/Debug mode)
    16 January 2005 -- added the adaptive delay based on the current width
    27 December 2004 -- initial version after experimenting all christmas long
  
===============================================================================
*/
#ifndef SMWRITER_SMD_H
#define SMWRITER_SMD_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"

class SMwriter_smd : public SMwriter
{
public:

  // additional mesh variables
  int nbits;

  // smwriter interface function implementations

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx);

  void close(bool close_file=true, bool update_header=true); // virtual base function must be reimplemented

  // smwriter_smd functions

  bool open(FILE* fd, int bits=16, int delay=-3);

  SMwriter_smd();
  ~SMwriter_smd();

private:
  int max_delay;
  int width_delay;

  void write_header();
  bool compress_triangle();
  bool compress_triangle_waiting();
  bool compress_triangle_traversal();
};

#endif
