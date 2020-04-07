/*
===============================================================================

  FILE:  smreader_smc_old.h
  
  CONTENTS:
  
    Reads a mesh from our compressed Streaming Mesh format (SMC) and provides
    access to it in form of a Streaming Mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    05 April 2005 -- slowly removing support for the old SMC reader and writer
    16 January 2004 -- added support for pre-existing bounding box info
    30 October 2003 -- switched to enums and bools in peter's office
    07 October 2003 -- initial version created the day of the California recall
  
===============================================================================
*/
#ifndef SMREADER_SMC_OLD_H
#define SMREADER_SMC_OLD_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_smc_old : public SMreader
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

  // smreader_smc_old functions

  bool open(FILE* file);

  SMreader_smc_old();
  ~SMreader_smc_old();

private:
  FILE* file;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  void read_header();
  int decompress_triangle();
};

#endif
