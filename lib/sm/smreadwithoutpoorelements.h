/*
===============================================================================

  FILE:  smreadwithoutpoorelements.h
  
  CONTENTS:
  
    A filter for reading a streaming mesh without those triangles whose ideal
    weight inverse mean ratio is above the cut_off (and also without vertices
    that have lost all their triangles).
 
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
#ifndef SMREAD_WITHOUT_POOR_ELEMENTS_H
#define SMREAD_WITHOUT_POOR_ELEMENTS_H

#include "smreaddereferenced.h"

#include <stdio.h>

class SMreadWithoutPoorElements : public SMreadDereferenced
{
public:

  // SMreadDereferenced interface function implementations

  virtual SMevent read_element();
  virtual SMevent read_event();

  virtual void close(bool close_file=true); // virtual base function must be reimplemented

  // SMreadWithoutPoorElements functions

  bool open(SMreader* smreader, float cut_off = 100.0f);

  SMreadWithoutPoorElements();
  ~SMreadWithoutPoorElements();

private:
  SMreader* smreader;
  int num_non_pre_order;
  float cut_off;

  int have_new, next_new;
  int new_vertices[3];
  int have_triangle;
  int have_finalized, next_finalized;
  int finalized_vertices[3];

  int read_input();
  int read_triangle();
};

#endif
