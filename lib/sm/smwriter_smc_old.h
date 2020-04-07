/*
===============================================================================

  FILE:  smwriter_smc_old.h
  
  CONTENTS:
  
    Writes a Streaming Mesh in a losslessly compressed binary format using a
    fairly "light-weight" compression scheme for the connectivity.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    05 April 2005 -- slowly removing support for the old SMC reader and writer
    3 January 2004 -- added support for pre-existing bounding box info
    28 October 2003 -- switched from Peter's hash to the STL hash_map
    15 September 2003 -- initial version created on the Monday after a tough
                         hiking weekend in Yosemite
  
===============================================================================
*/
#ifndef SMWRITER_SMC_OLD_H
#define SMWRITER_SMC_OLD_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"

class SMwriter_smc_old : public SMwriter
{
public:

  // smwriter interface function implementations

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx);

  void close(bool close_file=true, bool update_header=true); // virtual base function must be reimplemented

  // smwriter_smc_old functions

  bool open(FILE* fd, int bits=16);

  SMwriter_smc_old();
  ~SMwriter_smc_old();

private:
  FILE* file;
  void write_header();
};

#endif
