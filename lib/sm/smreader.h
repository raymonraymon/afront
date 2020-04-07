/*
===============================================================================

  FILE:  smreader.h
  
  CONTENTS:
  
    Streaming Mesh Reader
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003-07 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    28 May 2007 -- added constructor for initialization in derived classes
    23 May 2007 -- added the compute_bounding_box() function
    17 January 2004 -- added virtual destructor to shut up the g++ compiler
    30 October 2003 -- switched to enums and bools in peter's office
    02 October 2003 -- added functionality for read_element() 
    28 August 2003 -- initial version created on a hot Thursday at LLNL
  
===============================================================================
*/
#ifndef SMREADER_H
#define SMREADER_H

#include "smdefs.h"


class SMreader
{
public:
  // vertex variables
  int v_idx;
  float v_pos_f[3];

  // triangle variables
  int t_idx[3];
  bool t_final[3];

  // finalized variables
  int final_idx;

  // mesh variables

  int ncomments;
  char** comments;

  int nverts;
  int nfaces;

  int v_count;
  int f_count;

  float* bb_min_f;
  float* bb_max_f;

  bool post_order;

  // functions

  virtual bool compute_bounding_box(); // returns true if it leaves stream at end

  virtual SMevent read_element()=0;
  virtual SMevent read_event()=0;

  virtual void close(bool close_file=true);

  SMreader();
  virtual ~SMreader();
};

#endif
