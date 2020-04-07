/*
===============================================================================

  FILE:  smreader_jrs.h
  
  CONTENTS:
  
    Reads a mesh from an Jonathan R. Shewchuk's *.node and *.ele format and
    presents it as a (poor) streaming mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    24 May 2007 -- redefined more efficient compute_bounding_box() function
    24 August 2005 -- created after Henna decided not to come to Berkeley yet
  
===============================================================================
*/
#ifndef SMREADER_JRS_H
#define SMREADER_JRS_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_jrs : public SMreader
{
public:

  // smreader interface function implementations

  bool compute_bounding_box(); // reimplements virtual base function

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_jrs functions

  bool open(FILE* file_node, FILE* file_ele);

  SMreader_jrs();
  ~SMreader_jrs();

private:
  int ncoords;
  int skipped_lines;
  FILE* file_node;
  FILE* file_ele;
  char* line_node;
  char* line_ele;
};

#endif
