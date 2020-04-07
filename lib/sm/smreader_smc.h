/*
===============================================================================

  FILE:  smreader_smc.h
  
  CONTENTS:
  
    Reads a mesh from our improved compressed Streaming Mesh format (SMC) and provides
    access to it in form of a Streaming Mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 May 2005 -- fixed a Microsoft bug (floating-point in Release/Debug mode)
    05 April 2005 -- finally renamed from SMreader_sme to SMreader_smc
    21 March 2005 -- read_element() calls after EOF will always return SM_EOF 
    09 January 2005 -- temporary renamed (smf->sme), but will become the new smc
    07 January 2005 -- improved cache hits ratio (sme->smf)
    06 January 2005 -- further simplified and improved in speed (smd->sme)
    28 July 2004 -- improved in speed and rates using fewer tables (smc->smd)
    30 October 2003 -- switched to enums and bools in peter's office
    07 October 2003 -- initial version created the day of the California recall
  
===============================================================================
*/
#ifndef SMREADER_SMC_H
#define SMREADER_SMC_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_smc : public SMreader
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

  // smreader_smc functions

  bool open(FILE* file);

  SMreader_smc();
  ~SMreader_smc();

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
