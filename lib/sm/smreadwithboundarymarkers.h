/*
===============================================================================

  FILE:  smreadwithboundarymarkers.h
  
  CONTENTS:
  
    A filter for reading a streaming mesh with all boundary vertices marked.
    What the filter really checks for is whether a vertex is connected via
    the half-edges of *two* triangles to all its neighbours.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2007 martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    13 September 2007 -- after Germany beat Argentine 11:0 in Women Soccer 
  
===============================================================================
*/
#ifndef SMREAD_WITH_BOUNDARY_MARKERS_H
#define SMREAD_WITH_BOUNDARY_MARKERS_H

#include "smreaddereferenced.h"

#include <stdio.h>

class SMreadWithBoundaryMarkers : public SMreadDereferenced
{
public:
  // additional vertex variable
  bool v_boundary;

  // additional triangle variable
  bool t_boundary[3];

  // smreaddereferenced interface function implementations

  virtual SMevent read_element();
  virtual SMevent read_event();

  virtual void close(bool close_file=true); // virtual base function must be reimplemented

  // smreadwithboundarymarkers functions

  bool open(SMreader* smreader);

  SMreadWithBoundaryMarkers();
  ~SMreadWithBoundaryMarkers();

private:
  SMreader* smreader;
  int num_non_pre_order;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  int read_input();
  int read_triangle();
};

#endif
