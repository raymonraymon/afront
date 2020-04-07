/*
===============================================================================

  FILE:  smreadpreascompactpre.h
  
  CONTENTS:
  
    Reads a *pre-order* Streaming Mesh as a *compact* *pre-order* mesh. This
    requires reordering and reindexing of the vertices. It will fail if more
    than 'limit_buffer_size' triangles need to be stored. The 'preserve_order'
    flag specifies whether the algorithm is allowed to reorder triangles.

    If the algorithm must preserve the triangle order then the (triangle-)span
    of the stream dictates the maximal number of triangles it needs to buffer.
    This number will be prohibitively large for streaming meshes that were 
    created with some kind of depth first traversal.
    If the algorithm is allowed to reorder triangles, then the (vertex-)width
    of the stream roughly dictates the maximal number of triangles it needs to
    buffer.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 January 2005 -- major bug fix for behaviour after smreader had SM_EOF
     6 January 2004 -- setting the 'limit_buffer_size' may fail ungraceful
    25 December 2003 -- realized that there might be a bug for 'preserve_order == false'
    18 December 2003 -- refining the code after coming back to LLNL 
    10 October 2003 -- initial version created once SM compression worked
  
===============================================================================
*/
#ifndef SMREAD_PRE_AS_COMPACT_PRE_H
#define SMREAD_PRE_AS_COMPACT_PRE_H

#include "smreader.h"

class SMreadPreAsCompactPre : public SMreader
{
public:
  // additional triangle variables
  float* t_pos_f[3];

  // smreader interface function implementations

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // SMreadPreAsCompactPre functions

  bool open(SMreader* smreader, bool preserve_order=true, int limit_buffer_size=0);

  SMreadPreAsCompactPre();
  ~SMreadPreAsCompactPre();

private:
  SMreader* smreader;
  bool preserve_order;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  int read_preorder_input();
  int read_triangle();
};

#endif
