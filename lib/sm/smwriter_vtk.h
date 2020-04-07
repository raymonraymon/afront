/*
===============================================================================

  FILE:  smwriter_vtk.h
  
  CONTENTS:
  
    Writes a Streaming Triangle Mesh in the non-streaming VTK format.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    02 October 2007 -- created after depositing my first LLNL paycheck

===============================================================================
*/
#ifndef SMWRITER_VTK_H
#define SMWRITER_VTK_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"

class SMwriter_vtk : public SMwriter
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

  // smwriter_vtk functions

  bool open(FILE* file);

  SMwriter_vtk();
  ~SMwriter_vtk();

private:
  FILE* file;
  int vertex_buffer_alloc;
  float* vertex_buffer;
  int triangle_buffer_alloc;
  int* triangle_buffer;
};

#endif
