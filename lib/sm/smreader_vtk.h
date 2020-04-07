/*
===============================================================================

  FILE:  smreader_vtk.h
  
  CONTENTS:
  
    Reads a mesh from the VTK format provided the file contains a triangle mesh
    and presents it as a (poor) streaming mesh. Only the vertex positions and
    the triangle connectivity are extracted from the VTK file. Other fields are
    ignored.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003-07  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    06 September 2007 -- created at PanamaBay with acoustic "Everlong" live
  
===============================================================================
*/
#ifndef SMREADER_VTK_H
#define SMREADER_VTK_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_vtk : public SMreader
{
public:

  // smreader interface function implementations

  bool compute_bounding_box(); // reimplements virtual base function

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_vtk functions

  bool open(FILE* file);

  SMreader_vtk();
  ~SMreader_vtk();

private:
  int skipped_lines;
  FILE* file;
  char* line;
};

#endif
