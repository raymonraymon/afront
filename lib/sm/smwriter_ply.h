/*
===============================================================================

  FILE:  smwriter_ply.h
  
  CONTENTS:
  
    Writes a Streaming Mesh in standard PLY format. Since the output format PLY
    is not streaming all vertices and triangles have to be buffered. Hence, only
    small meshes can be converted like this.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    22 May 2007 -- created after Alex and Tamara invaded my Strada office

===============================================================================
*/
#ifndef SMWRITER_PLY_H
#define SMWRITER_PLY_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"

class SMwriter_ply : public SMwriter
{
public:

  // smwriter interface function implementations

  void set_nverts(int nverts);
  void set_nfaces(int nfaces);

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx);

  void close(bool close_file=true, bool update_header=true); // virtual base function must be reimplemented

  // smwriter_ply functions

  bool open(FILE* file);

  SMwriter_ply();
  ~SMwriter_ply();

private:
  void* out_ply;
  void write_header();
  int vertex_buffer_alloc;
  float* vertex_buffer;
  int triangle_buffer_alloc;
  int* triangle_buffer;
};

#endif
