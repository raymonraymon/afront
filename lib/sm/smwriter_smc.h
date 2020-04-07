/*
===============================================================================

  FILE:  smwriter_smc.h
  
  CONTENTS:
  
    Writes a Streaming Mesh in a improved compressed format using a stream-order
    preserving connectivity coding scheme. 
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2005  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    26 May 2005 -- fixed a Microsoft bug (floating-point in Release/Debug mode)
    05 April 2005 -- finally renamed from SMwriter_sme to SMwriter_smc
    04 April 2005 -- fixed mean SMC_END bug after SM paper was rejected ... again 
    09 January 2005 -- temporary renamed (smf->sme), but will become the new smc
    07 January 2005 -- improved cache hits ratio (sme->smf)
    05 January 2005 -- further simplified and improved in speed (smd->sme)
    27 July 2004 -- improved in speed and rates using fewer tables (smc->smd)
    3 January 2004 -- added support for pre-existing bounding box info
    28 October 2003 -- switched from Peter's hash to the STL hash_map
    15 September 2003 -- initial version created on the Monday after a tough
                         hiking weekend in Yosemite
  
===============================================================================
*/
#ifndef SMWRITER_SMC_H
#define SMWRITER_SMC_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"

class SMwriter_smc : public SMwriter
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

  // smwriter_smc functions

  bool open(FILE* fd, int bits=16);

  SMwriter_smc();
  ~SMwriter_smc();

private:
  void write_header();
};

#endif
