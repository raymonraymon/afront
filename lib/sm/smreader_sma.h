/*
===============================================================================

  FILE:  smreader_sma.h
  
  CONTENTS:
  
    Reads a mesh from our ASCII Streaming Mesh format (SMA) and provides
    access to it in form of a Streaming Mesh.
  
  PROGRAMMERS:
  
    martin isenburg@cs.unc.edu
  
  COPYRIGHT:
  
    copyright (C) 2003  martin isenburg@cs.unc.edu
    
    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
    15 January 2005 -- fixed valerio's bug (annoying output for empty lines) 
    07 January 2004 -- closing no longer resets comments, nverts, nfaces, bb_min
                       and bb_max. values are still there when re-opening again.
    19 December 2003 -- fixed bug: substract 1 from read finalized indices
    30 October 2003 -- switched to enums and bools in peter's office
    02 October 2003 -- added functionality for read_element() 
    28 August 2003 -- initial version created on a cold Thursday night
  
===============================================================================
*/
#ifndef SMREADER_SMA_H
#define SMREADER_SMA_H

#include "smreader.h"

#include <stdio.h>
#include "filewrapper.h"

class SMreader_sma : public SMreader
{
public:

  // smreader interface function implementations

  SMevent read_element();
  SMevent read_event();

  void close(bool close_file=true); // virtual base function must be reimplemented

  // smreader_sma functions

  bool open(FILE* fp);

  SMreader_sma();
  ~SMreader_sma();

private:
  FILE* file;
  int skipped_lines;
  char* line;
  int have_finalized, next_finalized;
  int finalized_vertices[3];
};

#endif
