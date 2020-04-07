/*
===============================================================================

  FILE:  smwriter_sma.h
  
  CONTENTS:
  
    Writes a Streaming Mesh in a simple ASCII format.
    
    These things we meant to do but have not gotten around to it yet: Optionally
    we would like the mesh indices to be relative or dynamic (pre-order only).
    When the mesh is written relative all indices would be expressed as the
    difference to the highest index. When the mesh is written dynamix then
    indices are re-used and possinble chance over time such that the maximal
    index used equals the cutwidth.
    
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    17 October 2007 -- changed way header is written so it can be "updated"
    02 October 2003 -- initial version created on the Thursday that Germany
                       beat Russia 7:1 in the Women Soccer Worlcup
  
===============================================================================
*/
#ifndef SMWRITER_SMA_H
#define SMWRITER_SMA_H

#include "smwriter.h"

#include <stdio.h>
#include "filewrapper.h"


class SMwriter_sma : public SMwriter
{
public:

  // smwriter interface function implementations

  void write_vertex(const float* v_pos_f);
  void write_triangle(const int* t_idx, const bool* t_final);
  void write_triangle(const int* t_idx);
  void write_finalized(int final_idx);

  void close(bool close_file=true, bool update_header=true); // virtual base function must be reimplemented

  // smwriter_sma functions

  bool open(FILE* file);

  SMwriter_sma();
  ~SMwriter_sma();

private:
  FILE* file;
  void write_header();
  void update_header();
};

#endif
