/*
===============================================================================

  FILE:  smreader_smd.h
  
  CONTENTS:
  
    Reads a compressed Streaming Mesh using a more involved compression scheme
    that utilized a small buffer from which it chooses its triangles. There is
    a user-defined constraint on the maximum that a triangle can be delayed in
    the buffer.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 May 2005 -- fixed a Microsoft bug (floating-point in Release/Debug mode)
    21 March 2005 -- read_element() calls after EOF will always return SM_EOF 
    04 January 2005 -- created after returning on a red-eye from Calgary
  
===============================================================================
*/
#ifndef SMREADER_SMD_H
#define SMREADER_SMD_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_smd : public SMreader
{
public:
  // additional triangle variables
  float* t_pos_f[3];

  // additional mesh variables
  int nbits;

  // smreader interface function implementations

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_smx functions

  bool open(FILE* file);

  SMreader_smd();
  ~SMreader_smd();

private:
  FILE* file;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  void read_header();
  bool decompress_triangle();
  bool decompress_triangle_waiting();
  bool decompress_triangle_traversal();
};

#endif
