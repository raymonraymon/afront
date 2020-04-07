/*
===============================================================================

  FILE:  smreader_ply.h
  
  CONTENTS:
  
    Reads a vertices and triangles of a mesh in Stanford's non-streaming PLY
    format and provides access to it in form of a Streaming Mesh of maximal
    width.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    24 May 2007 -- redefined more efficient compute_bounding_box() function
    17 January 2004 -- v_idx was not being set correctlty 
    6 January 2004 -- added support to compute the bounding box
    22 December 2003 -- initial version created once christmas plans were off
  
===============================================================================
*/
#ifndef SMREADER_PLY_H
#define SMREADER_PLY_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_ply : public SMreader
{
public:

  // smreader interface function implementations

  bool compute_bounding_box(); // reimplements virtual base function

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_ply functions

  bool open(FILE* file);

  SMreader_ply();
  ~SMreader_ply();

private:
  void* in_ply;
  int velem;
  int felem;
};

#endif
